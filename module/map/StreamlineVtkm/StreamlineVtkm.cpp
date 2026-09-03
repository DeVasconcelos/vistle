#include <viskores/filter/flow/Streamline.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/VectorAnalysis.h>
#include <viskores/cont/EnvironmentTracker.h>
#include <viskores/thirdparty/diy/diy.h>

#include <vistle/util/enum.h>
#include <vistle/vtkm/convert.h>
#include <vistle/vtkm/vtkm_module_utils.h>

#include "worklet/GenerateSeeds.h"
#include "StreamlineVtkm.h"

using namespace vistle;

MODULE_MAIN(StreamlineVtkm)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(IntegrationMethod, (RK4)(Euler))
DEFINE_ENUM_WITH_STRING_CONVERSIONS(StartStyle, (Line)(Plane))

StreamlineVtkm::StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: Module(name, moduleID, comm), m_numPorts(3), m_mappedDataHandling(MappedDataHandling::Require)
{
    setReducePolicy(message::ReducePolicy::PerTimestep);

    const MPI_Comm mpiComm = comm;
    const viskoresdiy::mpi::DIY_MPI_Comm diyComm{mpiComm};
    viskoresdiy::mpi::communicator viskoresComm(diyComm, false);
    viskoresdiy::mpi::communicator viskoresDup;
    viskoresDup.duplicate(viskoresComm);

    viskores::cont::EnvironmentTracker::SetCommunicator(viskoresComm);

    assert(m_numPorts > 0);
    bool dataInput =
        m_mappedDataHandling != MappedDataHandling::Discard && m_mappedDataHandling != MappedDataHandling::Generate;
    bool dataOutput = m_mappedDataHandling != MappedDataHandling::Discard;

    for (int i = 0; i < m_numPorts; ++i) {
        std::string in("data_in");
        std::string out("data_out");
        if (i > 0) {
            in += std::to_string(i);
            out += std::to_string(i);
        }
        m_inputPorts.push_back(createInputPort(in, dataInput ? "input grid with mapped data" : "input grid"));
        m_outputPorts.push_back(createOutputPort(out, dataOutput ? "output grid with mapped data" : "output grid"));
        linkPorts(m_inputPorts[i], m_outputPorts[i]);
        if (i > 0) {
            setPortOptional(m_inputPorts[i], true);
        }
    }

    m_printObjectInfo =
        addIntParameter("_print_object_info", "print information on generated data objects for debug purposes", false,
                        Parameter::Boolean);
    // ------------------------------------------------------------------------------

    setCurrentParameterGroup("Seed Points");
    const Integer max_no_seeds = 300;
    m_numberOfSeeds = addIntParameter("number_of_seeds", "number of seed points", 2);
    setParameterRange(m_numberOfSeeds, (Integer)1, max_no_seeds);
    m_startStyle =
        addIntParameter("start_style", "initial particle position configuration", StartStyle::Line, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(m_startStyle, StartStyle, );
    m_startPoint1 = addVectorParameter("startpoint1", "1st initial point", ParamVector(0, 0.2, 0));
    m_startPoint2 = addVectorParameter("startpoint2", "2nd initial point", ParamVector(1, 0, 0));
    m_direction = addVectorParameter("direction", "tracing direction", ParamVector(0, 0, 1));
    m_maxNumberOfSeeds =
        addIntParameter("max_no_seeds", "maximum number of seeds (for parameter/slider limits)", max_no_seeds);
    setParameterRange(m_maxNumberOfSeeds, (Integer)2, (Integer)1000000);

    setCurrentParameterGroup("Stop Conditions");
    m_numberOfSteps = addIntParameter("steps_max", "maximum number of integration steps", 100);

    setCurrentParameterGroup("Step Length Control");
    m_integrationMethod =
        addIntParameter("integration_method", "integration method", IntegrationMethod::RK4, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(m_integrationMethod, IntegrationMethod, );
    m_stepSize = addFloatParameter("step_size", "integration step size", 0.1f);
}

std::string StreamlineVtkm::getFieldName(int index, bool output) const
{
    std::string name = "data_at_port_" + std::to_string(index);
    if (index == 0 && output)
        name += "_out";
    return name;
}

bool StreamlineVtkm::prepare()
{
    if (!m_inputPorts[0]->isConnected()) {
        if (rank() == 0)
            sendError("No input connected to %s", m_inputPorts[0]->getName().c_str());
        return false;
    }

    for (int i = 0; i < m_numPorts; ++i) {
        if (!m_inputPorts[i]->isConnected() && m_outputPorts[i]->isConnected()) {
            if (rank() == 0)
                sendError("Output port " + m_outputPorts[i]->getName() +
                          " is connected, but corresponding input port " + m_inputPorts[i]->getName() + " is not");
            return false;
        }
    }

    return Module::prepare();
}

ModuleStatusPtr StreamlineVtkm::readInPorts(const std::shared_ptr<BlockTask> &task, Object::const_ptr &grid,
                                            std::vector<DataBase::const_ptr> &fields) const
{
    for (int i = 0; i < m_numPorts; ++i) {
        if (!m_inputPorts[i]->isConnected()) {
            fields.push_back(nullptr);
            continue;
        }

        auto container = task->accept<Object>(m_inputPorts[i]);
        auto split = splitContainerObject(container);
        auto geometry = split.geometry;
        auto data = split.mapped;

        // make sure there is data on the input port if the corresponding output port is connected
        if (!geometry && !data)
            return Error("No data on input port " + m_inputPorts[i]->getName() + ", even though it is connected");

        fields.push_back(data);

        // make sure all data fields are defined on the same grid
        if (grid) {
            if (geometry && geometry->getHandle() != grid->getHandle()) {
                return Error("The grid on " + m_inputPorts[i]->getName() +
                             " does not match the grid on the other input ports!");
            }
        } else {
            grid = geometry;
        }
    }

    if (!grid)
        return Error("Could not find a valid input grid on any input port");

    return Success();
}

bool StreamlineVtkm::changeParameter(const Parameter *param)
{
    if (param == m_startStyle) {
        setParameterReadOnly(m_direction, m_startStyle->getValue() != Plane);
    } else if (param == m_maxNumberOfSeeds) {
        setParameterRange(m_numberOfSeeds, (Integer)1, m_maxNumberOfSeeds->getValue());
    }

    return Module::changeParameter(param);
}

ModuleStatusPtr StreamlineVtkm::prepareInputGrid(InputData &input) const
{
    return vtkmSetGrid(input.viskoresDataset, input.vistleGrid);
}

ModuleStatusPtr StreamlineVtkm::prepareInputField(const Port *port, InputData &input, int index) const
{
    auto field = input.fields[index];

    if (index == 0) {
        if (auto in = Vec<Scalar, 3>::as(field); !in) {
            return Error("Error: Input field at port " + port->getName() + " must be a 3D vector field!");
        }
    }

    return vtkmAddField(input.viskoresDataset, field, getFieldName(index));
}

bool StreamlineVtkm::reduce(int timestep)
{
    auto filter = setUpFilter();
    filter->SetActiveField(getFieldName(0));

    filter->SetOutputFieldName(getFieldName(0, true));
    filter->SetFieldsToPass("", viskores::cont::Field::Association::Any, viskores::filter::FieldSelection::Mode::All);

    viskores::cont::PartitionedDataSet output;
    if (!this->tryToExecuteFilter(*filter, m_globalData.partitionedDataset, output))
        return true;

    std::lock_guard<std::mutex> lock(m_globalData.mutex);
    m_globalData.outputGrids.clear();
    m_globalData.outputFields.assign(m_numPorts, {});

    for (viskores::Id partitionIndex = 0; partitionIndex < output.GetNumberOfPartitions(); ++partitionIndex) {
        const auto &dataset = output.GetPartition(partitionIndex);
        auto outputGrid = vtkmGetGeometry(dataset);
        if (!outputGrid) {
            sendError("Could not convert StreamlineVtkm output geometry to a Vistle object.");
            continue;
        }

        updateMeta(outputGrid);
        m_globalData.outputGrids.push_back(outputGrid);

        for (int port = 0; port < m_numPorts; ++port) {
            if (!m_outputPorts[port]->isConnected())
                continue;

            std::string fieldName = getFieldName(port);
            if (port == 0 && dataset.HasField(getFieldName(0, true)))
                fieldName = getFieldName(0, true);

            auto field = vtkmGetField(dataset, fieldName, DataBase::Unspecified, false);
            if (field) {
                updateMeta(field);
                field->setGrid(outputGrid);
                m_globalData.outputFields[port].push_back(field);
                addObject(m_outputPorts[port], field);
            } else {
                m_globalData.outputFields[port].push_back(nullptr);
                addObject(m_outputPorts[port], outputGrid);
            }
        }
    }

    return true;
}

bool StreamlineVtkm::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    InputData input;
    std::vector<std::vector<vistle::Object::const_ptr>> grid;
    std::vector<std::vector<std::vector<vistle::DataBase::const_ptr>>> data_in;
    std::vector<viskores::cont::DataSet> viskoresDatasets;

    auto status = readInPorts(task, input.vistleGrid, input.fields);
    if (!checkAndNotify(status))
        return true;

    assert(m_outputPorts.size() == input.fields.size());

    auto timestep = input.vistleGrid->getTimestep();
    if (timestep < 0 && input.fields[0])
        timestep = input.fields[0]->getTimestep();
    if (timestep < 0)
        timestep = -1;

    const auto timestepIndex = static_cast<std::size_t>(timestep + 1);
    if (grid.size() <= timestepIndex) {
        grid.resize(timestepIndex + 1);
        data_in.resize(m_numPorts);
        for (auto &portData: data_in)
            portData.resize(timestepIndex + 1);
    } else {
        for (auto &portData: data_in) {
            if (portData.size() <= timestepIndex)
                portData.resize(timestepIndex + 1);
        }
    }

    grid[timestepIndex].push_back(input.vistleGrid);
    viskoresDatasets.emplace_back();
    status = vtkmSetGrid(viskoresDatasets.back(), input.vistleGrid);
    if (!checkAndNotify(status))
        return true;

    for (std::size_t i = 0; i < input.fields.size(); ++i) {
        data_in[i][timestepIndex].push_back(input.fields[i]);
        auto field = data_in[i][timestepIndex].back();
        if (field) {
            status = vtkmAddField(viskoresDatasets.back(), field, getFieldName(i));
            if (!checkAndNotify(status))
                return true;
        }
    }

    std::lock_guard<std::mutex> lock(m_globalData.mutex);
    m_globalData.partitionedDataset.AppendPartitions(viskoresDatasets);

    return true;
#if 0
    InputData input;
    OutputData output;

    bool printInfo = m_printObjectInfo->getValue() != 0;

    // read in data from the input ports...
    auto status = readInPorts(task, input.vistleGrid, input.fields);
    if (!checkAndNotify(status))
        return true;

    assert(m_outputPorts.size() == input.fields.size());

    // ... transform the input grid (and fields) into a Viskores dataset ...
    status = prepareInputGrid(input);
    if (!checkAndNotify(status))
        return true;

    for (std::size_t i = 0; i < input.fields.size(); ++i) {
        if (i > 0 && !m_outputPorts[i]->isConnected())
            continue;

        if (m_mappedDataHandling == MappedDataHandling::Require) {
            if (!input.fields[i]) {
                sendError("Cannot continue: No mapped data on input port " + m_inputPorts[i]->getName() +
                          ", even though it is required by the module!");
                return true;
            }
        }
        if (input.fields[i]) {
            status = prepareInputField(m_inputPorts[i], input, i);
            if (!checkAndNotify(status))
                return true;
        }
    }

    // ... run filter on the active field ...
    bool useInputData =
        m_mappedDataHandling != MappedDataHandling::Discard && m_mappedDataHandling != MappedDataHandling::Generate;
    auto activeField = useInputData ? getFieldName(0) : "";

    if (printInfo) {
        std::stringstream str;
        str << "<pre>Input ";
        input.viskoresDataset.PrintSummary(str);
        str << "</pre>" << std::endl;
        std::cout << str.str() << std::endl;
    }

    if (m_mappedDataHandling != MappedDataHandling::Require || input.viskoresDataset.HasField(activeField)) {
        if (auto filter = setUpFilter()) {
            if (input.viskoresDataset.HasField(activeField))
                filter->SetActiveField(activeField);

            /*
                By default, Viskores names output fields the same as input fields which causes problems
                if the input mapping is different from the output mapping, i.e., when converting
                a point field to a cell field or vice versa. To avoid having a point and a
                cell field of the same name in the resulting dataset, which leads to conflicts, e.g.,
                when calling Viskores's GetField() method, we rename the output field here.
            */
            filter->SetOutputFieldName(getFieldName(0, true));
            filter->SetFieldsToPass("", viskores::cont::Field::Association::Any,
                                    viskores::filter::FieldSelection::Mode::All);

            if (printInfo) {
                std::stringstream str;
                str << "Filter: " << typeid(decltype(*filter)).name() << std::endl;
                std::cout << str.str() << std::endl;
            }

            if (!this->tryToExecuteFilter(*filter, input.viskoresDataset, output.viskoresDataset))
                return true;

            if (printInfo) {
                std::stringstream str;
                str << "<pre>Output ";
                output.viskoresDataset.PrintSummary(str);
                str << "</pre>" << std::endl;
                std::cout << str.str() << std::endl;
            }
        } else {
            output.viskoresDataset = input.viskoresDataset;

            if (printInfo) {
                std::stringstream str;
                str << "<pre>Output ";
                output.viskoresDataset.PrintSummary(str);
                str << "</pre>" << std::endl;
                std::cout << str.str() << std::endl;
            }
        }
    }

    // ... transform filter output, i.e., grid and data fields, to Vistle objects ...
    output.vistleGrid = prepareOutputGrid(input, output);
    if (!output.vistleGrid)
        return true;

    output.fields.resize(input.fields.size(), nullptr);
    for (std::size_t i = 0; i < input.fields.size(); ++i) {
        if (!m_outputPorts[i]->isConnected())
            continue;

        if (m_mappedDataHandling != MappedDataHandling::Use || input.fields[i]) {
            std::string outputFieldName = getFieldName(i);
            if (i == 0 && output.viskoresDataset.HasField(getFieldName(i, true))) {
                // if filter has created a dedicated output field, use it
                outputFieldName = getFieldName(i, true);
            }
            output.fields[i] = prepareOutputField(input, output, i, outputFieldName);
        }

        // ... and write the result to the output ports
        if (output.fields[i] || m_mappedDataHandling == MappedDataHandling::Generate) {
            task->addObject(m_outputPorts[i], output.fields[i]);
        } else if (output.vistleGrid) {
            task->addObject(m_outputPorts[i], output.vistleGrid);
        }
    }

    return true;
#endif
}

viskores::cont::ArrayHandle<viskores::Particle> StreamlineVtkm::createSeedArray() const
{
    viskores::Vec3f startpoint1{m_startPoint1->getValue()[0], m_startPoint1->getValue()[1],
                                m_startPoint1->getValue()[2]};
    viskores::Vec3f startpoint2{m_startPoint2->getValue()[0], m_startPoint2->getValue()[1],
                                m_startPoint2->getValue()[2]};

    if (m_startStyle->getValue() == StartStyle::Line) {
        return generateSeedsOnLine(m_numberOfSeeds->getValue(), startpoint1, startpoint2);

    } else {
        viskores::Vec3f direction{m_direction->getValue()[0], m_direction->getValue()[1], m_direction->getValue()[2]};
        return generateSeedsOnPlane(m_numberOfSeeds->getValue(), startpoint1, startpoint2, direction);
    }
}

std::unique_ptr<viskores::filter::Filter> StreamlineVtkm::setUpFilter() const
{
    auto filter = std::make_unique<viskores::filter::flow::Streamline>();

    filter->SetStepSize(m_stepSize->getValue());
    filter->SetNumberOfSteps(m_numberOfSteps->getValue());

    auto seedArray = createSeedArray();
    filter->SetSeeds(seedArray);

    if (m_integrationMethod->getValue() == IntegrationMethod::Euler)
        filter->SetSolverEuler();
    else
        filter->SetSolverRK4();

    return filter;
}

Object::const_ptr StreamlineVtkm::prepareOutputGrid(const InputData &input, OutputData &output) const
{
    auto outputGrid = vtkmGetGeometry(output.viskoresDataset);
    if (outputGrid) {
        updateMeta(outputGrid);
        outputGrid->copyAttributes(input.vistleGrid);
    }

    return outputGrid;
}

DataBase::ptr StreamlineVtkm::prepareOutputField(const InputData &input, OutputData &output, int index,
                                                 const std::string &fieldName) const
{
    auto probe = viskores::filter::resampling::Probe();
    probe.SetGeometry(output.viskoresDataset);
    probe.SetOutputFieldName(fieldName);

    viskores::cont::DataSet probeOutput;
    if (!this->tryToExecuteFilter(probe, input.viskoresDataset, probeOutput)) {
        sendError("An error occurred while probing the filter output field " + fieldName + " to the output grid.");
        return nullptr;
    }

    // --------------------------------------------------------------------------------

    if (auto mapped = vtkmGetField(probeOutput, fieldName)) {
        std::cerr << "mapped data: " << *mapped << std::endl;
        updateMeta(mapped);

        // the mapping of the output field might differ from the one of the input field,
        // so lets temporarily store it and set it again after copying the attributes
        auto mapping = mapped->mapping();
        mapped->copyAttributes(input.fields[index]);
        mapped->setMapping(mapping);

        if (output.vistleGrid)
            mapped->setGrid(output.vistleGrid);

        return mapped;

    } else {
        sendError("An error occurred while transforming the filter output field " + fieldName + " to a Vistle object.");
    }

    return nullptr;
}

bool StreamlineVtkm::tryToExecuteFilter(viskores::filter::Filter &filter, const viskores::cont::DataSet &inputDataset,
                                        viskores::cont::DataSet &outputDataset) const
{
    return vistle::vtkm::tryToExecuteFilter(*this, filter, inputDataset, outputDataset);
}

bool StreamlineVtkm::tryToExecuteFilter(viskores::filter::Filter &filter,
                                        const viskores::cont::PartitionedDataSet &inputDataset,
                                        viskores::cont::PartitionedDataSet &outputDataset) const
{
    return vistle::vtkm::tryToExecuteFilter(*this, filter, inputDataset, outputDataset);
}

bool StreamlineVtkm::checkAndNotify(const ModuleStatusPtr &status) const
{
    return vistle::vtkm::checkAndNotify(*this, status);
}
