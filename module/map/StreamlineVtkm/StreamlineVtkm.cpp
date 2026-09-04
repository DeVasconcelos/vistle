#include <viskores/filter/flow/Streamline.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/VectorAnalysis.h>
#include <viskores/cont/EnvironmentTracker.h>
#include <viskores/cont/MergePartitionedDataSet.h>
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

void StreamlineVtkm::GlobalData::clear()
{
    partitionedDatasets.clear();
    inputGrids.clear();
    inputFields.clear();
}

bool StreamlineVtkm::GlobalData::isEmpty() const
{
    return partitionedDatasets.empty() && inputGrids.empty() && inputFields.empty();
}

void StreamlineVtkm::GlobalData::resize(std::size_t newSize, std::size_t numFields)
{
    if (inputFields.size() < numFields)
        inputFields.resize(numFields);

    if (partitionedDatasets.size() < newSize) {
        partitionedDatasets.resize(newSize);
        inputGrids.resize(newSize);
        for (auto &field: inputFields)
            field.resize(newSize);
    }
}

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

    std::lock_guard<std::mutex> lock(m_globalData.mutex);
    if (!m_globalData.isEmpty())
        m_globalData.clear();

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

bool StreamlineVtkm::reduce(int timestep)
{
    const auto idx = static_cast<std::size_t>(timestep + 1);

    viskores::cont::PartitionedDataSet inputPartitionedDataset;
    std::vector<vistle::Object::const_ptr> inputGrids;
    std::vector<std::vector<vistle::DataBase::const_ptr>> inputFieldsForPorts; // [port][partition]
    {
        std::lock_guard<std::mutex> lock(m_globalData.mutex);
        // reduce(-1) is always called once in addition to reduce() for every actual timestep (see
        // Module::reduceWrapper), so it is normal for there to be no data left to process for it
        if (idx >= m_globalData.partitionedDatasets.size())
            return true;
        inputPartitionedDataset = m_globalData.partitionedDatasets[idx];
        if (inputPartitionedDataset.GetNumberOfPartitions() == 0)
            return true;

        if (idx < m_globalData.inputGrids.size())
            inputGrids = m_globalData.inputGrids[idx];

        inputFieldsForPorts.resize(m_globalData.inputFields.size());
        for (std::size_t port = 0; port < m_globalData.inputFields.size(); ++port) {
            if (idx < m_globalData.inputFields[port].size())
                inputFieldsForPorts[port] = m_globalData.inputFields[port][idx];
        }
    }

    auto filter = setUpFilter();
    filter->SetActiveField(getFieldName(0));

    filter->SetOutputFieldName(getFieldName(0, true));
    filter->SetFieldsToPass("", viskores::cont::Field::Association::Any, viskores::filter::FieldSelection::Mode::All);

    viskores::cont::PartitionedDataSet output;
    if (!this->tryToExecuteFilter(*filter, inputPartitionedDataset, output))
        return true;

    // the Streamline filter can drop input partitions that never received a particle and does not
    // otherwise report which input partition an output partition came from, so the output partitions
    // cannot be positionally matched to the input ones; since streamlines can also cross block
    // boundaries, merge all blocks of this timestep into a single dataset to probe fields from
    auto mergedInputDataset = viskores::cont::MergePartitionedDataSet(inputPartitionedDataset);

    // like Tracer, tag the generated grid/fields with the timestep they belong to
    Meta meta;
    meta.setNumBlocks(size());
    meta.setBlock(rank());
    meta.setNumTimesteps(numTimesteps() > 0 ? numTimesteps() : -1);
    meta.setTimeStep(numTimesteps() > 0 ? timestep : -1);

    for (viskores::Id partitionIndex = 0; partitionIndex < output.GetNumberOfPartitions(); ++partitionIndex) {
        const auto &dataset = output.GetPartition(partitionIndex);
        auto outputGrid = vtkmGetGeometry(dataset);
        if (!outputGrid) {
            sendError("Could not convert StreamlineVtkm output geometry to a Vistle object.");
            continue;
        }

        // setMeta() must come before updateMeta(), since it resets creator/generation to defaults,
        // which updateMeta() then fills in correctly; copyAttributes() only affects the attribute list
        outputGrid->setMeta(meta);
        updateMeta(outputGrid);
        // attributes (e.g. species name, color map range) are the same on every block of a given
        // timestep, so any input grid can be used as the source, regardless of which output partition
        // this is (output partitions cannot be positionally matched to input ones, see above)
        if (!inputGrids.empty()) {
            if (auto &vistleGrid = inputGrids.front())
                outputGrid->copyAttributes(vistleGrid);
        }

        for (int port = 0; port < m_numPorts; ++port) {
            if (!m_outputPorts[port]->isConnected())
                continue;

            DataBase::ptr field;
            if (port == 0 && dataset.HasField(getFieldName(0, true))) {
                // the field that was integrated is already attached to the lines by the Streamline filter
                field = vtkmGetField(dataset, getFieldName(0, true), DataBase::Unspecified, false);
            } else {
                // the Streamline filter does not map the other input fields onto the generated lines,
                // so we have to probe them from the merged input dataset
                auto probe = viskores::filter::resampling::Probe();
                probe.SetGeometry(dataset);
                probe.SetOutputFieldName(getFieldName(port));

                viskores::cont::DataSet probeOutput;
                if (this->tryToExecuteFilter(probe, mergedInputDataset, probeOutput)) {
                    field = vtkmGetField(probeOutput, getFieldName(port));
                }
            }

            if (field) {
                // setMeta() must come before updateMeta(), for the same reason as for outputGrid above
                field->setMeta(meta);
                updateMeta(field);
                vistle::DataBase::const_ptr origField;
                if (static_cast<std::size_t>(port) < inputFieldsForPorts.size() && !inputFieldsForPorts[port].empty())
                    origField = inputFieldsForPorts[port].front();
                if (origField) {
                    auto mapping = field->mapping();
                    field->copyAttributes(origField);
                    field->setMapping(mapping);
                }
                field->setGrid(outputGrid);
                addObject(m_outputPorts[port], field);
            } else {
                addObject(m_outputPorts[port], outputGrid);
            }
        }
    }

    return true;
}

bool StreamlineVtkm::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    vistle::Object::const_ptr vistleGrid;
    std::vector<vistle::DataBase::const_ptr> fields;

    viskores::cont::DataSet viskoresDataset;

    auto status = readInPorts(task, vistleGrid, fields);
    if (!checkAndNotify(status))
        return true;

    assert(m_outputPorts.size() == fields.size());

    auto timestep = vistleGrid->getTimestep();

    if (fields[0]) {
        if (timestep != fields[0]->getTimestep()) {
            sendError("timestep mismatch: grid = %d, field = %d", timestep, fields[0]->getTimestep());
            return true;
        }

        if (auto in = Vec<Scalar, 3>::as(fields[0]); !in) {
            sendError("Error: Input field at port " + m_inputPorts[0]->getName() + " must be a 3D vector field!");
            return true;
        }
    }

    if (timestep < 0)
        timestep = -1;

    status = vtkmSetGrid(viskoresDataset, vistleGrid);
    if (!checkAndNotify(status))
        return true;

    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (fields[i]) {
            status = vtkmAddField(viskoresDataset, fields[i], getFieldName(i));
            if (!checkAndNotify(status))
                return true;
        }
    }

    std::lock_guard<std::mutex> lock(m_globalData.mutex);
    auto numSteps = static_cast<std::size_t>(timestep + 1);
    m_globalData.resize(numSteps + 1, fields.size());

    m_globalData.partitionedDatasets[numSteps].AppendPartitions({viskoresDataset});
    m_globalData.inputGrids[numSteps].push_back(vistleGrid);
    for (std::size_t i = 0; i < fields.size(); ++i)
        m_globalData.inputFields[i][numSteps].push_back(fields[i]);

    return true;
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
