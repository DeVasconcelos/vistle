#include <viskores/cont/Error.h>
#include <viskores/cont/ErrorBadAllocation.h>
#include <viskores/cont/ErrorBadDevice.h>
#include <viskores/cont/ErrorBadType.h>
#include <viskores/cont/ErrorBadValue.h>
#include <viskores/cont/ErrorExecution.h>
#include <viskores/cont/ErrorFilterExecution.h>
#include <viskores/cont/ErrorInternal.h>

#include <viskores/filter/flow/Streamline.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/VectorAnalysis.h>

#include <vistle/util/enum.h>
#include <vistle/vtkm/convert.h>

#include "worklet/GenerateSeeds.h"
#include "StreamlineVtkm.h"

using namespace vistle;

MODULE_MAIN(StreamlineVtkm)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(IntegrationMethod, (RK4)(Euler))
DEFINE_ENUM_WITH_STRING_CONVERSIONS(StartStyle, (Line)(Plane))

StreamlineVtkm::StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: Module(name, moduleID, comm), m_numPorts(3), m_mappedDataHandling(MappedDataHandling::Require)
{
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

std::string StreamlineVtkm::getFieldName(int i, bool output) const
{
    std::string name = "data_at_port_" + std::to_string(i);
    if (i == 0 && output)
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
            std::stringstream msg;
            msg << "Output port " << m_outputPorts[i]->getName() << " is connected, but corresponding input port "
                << m_inputPorts[i]->getName() << " is not";
            if (rank() == 0)
                sendError(msg.str());
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
        if (!geometry && !data) {
            std::stringstream msg;
            msg << "No data on input port " << m_inputPorts[i]->getName() << ", even though it is connected";
            return Error(msg.str());
        }

        fields.push_back(data);

        // make sure all data fields are defined on the same grid
        if (grid) {
            if (geometry && geometry->getHandle() != grid->getHandle()) {
                std::stringstream msg;
                msg << "The grid on " << m_inputPorts[i]->getName()
                    << " does not match the grid on the other input ports!";
                return Error(msg.str());
            }
        } else {
            grid = geometry;
        }
    }

    if (!grid) {
        std::ostringstream msg;
        msg << "Could not find a valid input grid on any input port";
        return Error(msg.str());
    }

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
    /*
        TODO: 
        As the module is now, it cannot handle partitioned input grids as the Streamline
        filter will not know that it has to exchange particles with other partitions.
        While it seems like the filter does support this, given a viskores::cont::Par-
        titionedDataSet (see viskores/viskores/filter/flow/internal/ParticleExchanger.h),
        as far as I understand, for us, each block creates their own unpartitioned dataset.
        The filter either does so serially or using MPI (if VISKORES_ENABLE_MPI is ON, which,
        in our case, it never is).

        So, we could either:
        - implement a custom filter that handles the particle exchange between partitions 
          (similar to the one Viskores has, but for our setup). This would also require digging
          into how MPI communication is handled in Viskores, especially when GPUs are involved
        - when converting a partitioned block to the Viskores data format, use viskores::cont::
          PartitionedDataSet, although this would affect a lot of code, and in the end, still
          be serial, as VISKORES_ENABLE_MPI is OFF
        - have the master rank create a viskores::cont::PartitionedDataSet here (although again
          this would be serial, as VISKORES_ENABLE_MPI is OFF)
    */
    if (input.vistleGrid->getNumBlocks() != 1)
        return Error("StreamlineVtkm: Partitioned input grids are not supported yet!");

    return vtkmSetGrid(input.viskoresDataset, input.vistleGrid);
}

ModuleStatusPtr StreamlineVtkm::prepareInputField(const Port *port, InputData &input, int index) const
{
    auto field = input.fields[index];

    if (port->getName() == "data_in") {
        if (auto in = Vec<Scalar, 3>::as(field); !in) {
            return Error("Error: Input field at port " + port->getName() + " must be a 3D vector field!");
        }
    }

    return vtkmAddField(input.viskoresDataset, field, getFieldName(index));
}

bool StreamlineVtkm::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    InputData input;
    OutputData output;

    bool printInfo = m_printObjectInfo->getValue() != 0;

    // read in data from the input ports...
    auto status = readInPorts(task, input.vistleGrid, input.fields);
    if (!isValid(status))
        return true;

    assert(m_outputPorts.size() == input.fields.size());

    // ... transform the input grid (and fields) into a Viskores dataset ...
    status = prepareInputGrid(input);
    if (!isValid(status))
        return true;

    for (std::size_t i = 0; i < input.fields.size(); ++i) {
        if (i > 0 && !m_outputPorts[i]->isConnected())
            continue;

        if (m_mappedDataHandling == MappedDataHandling::Require) {
            if (!input.fields[i]) {
                std::stringstream msg;
                msg << "No mapped data on input port " << m_inputPorts[i]->getName();
                status = Error(msg.str());
                if (!isValid(status))
                    return true;
            }
        }
        if (input.fields[i]) {
            status = prepareInputField(m_inputPorts[i], input, i);
            if (!isValid(status))
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
        auto msg = str.str();
        std::cout << msg << std::endl;
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
                auto msg = str.str();
                std::cout << msg << std::endl;
            }

            if (!tryToExecuteFilter(filter, input.viskoresDataset, output.viskoresDataset))
                return true;

            if (printInfo) {
                std::stringstream str;
                str << "<pre>Output ";
                output.viskoresDataset.PrintSummary(str);
                str << "</pre>" << std::endl;
                auto msg = str.str();
                std::cout << msg << std::endl;
            }
        } else {
            output.viskoresDataset = input.viskoresDataset;

            if (printInfo) {
                std::stringstream str;
                str << "<pre>Output ";
                output.viskoresDataset.PrintSummary(str);
                str << "</pre>" << std::endl;
                auto msg = str.str();
                std::cout << msg << std::endl;
            }
        }
    }

    // ... transform filter output, i.e., grid and data fields, to Vistle objects ...
    status = prepareOutputGrid(input, output);
    if (!isValid(status))
        return true;

    output.fields.resize(input.fields.size(), nullptr);
    for (std::size_t i = 0; i < input.fields.size(); ++i) {
        if (!m_outputPorts[i]->isConnected())
            continue;

        if (m_mappedDataHandling != MappedDataHandling::Use || input.fields[i]) {
            status = prepareOutputField(input, output, i);
            if (!isValid(status))
                return true;
        }

        // ... and write the result to the output ports
        if (output.fields[i] || m_mappedDataHandling == MappedDataHandling::Generate) {
            task->addObject(m_outputPorts[i], output.fields[i]);
        } else if (output.vistleGrid) {
            task->addObject(m_outputPorts[i], output.vistleGrid);
        }
    }

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

ModuleStatusPtr StreamlineVtkm::prepareOutputGrid(const InputData &input, OutputData &output) const
{
    output.vistleGrid = vtkmGetGeometry(output.viskoresDataset);
    if (output.vistleGrid) {
        updateMeta(output.vistleGrid);
        output.vistleGrid->copyAttributes(input.vistleGrid);
    } else {
        return Error("An error occurred while transforming the filter output grid to a Vistle object.");
    }

    return Success();
}

ModuleStatusPtr StreamlineVtkm::prepareOutputField(const InputData &input, OutputData &output, int index) const
{
    std::string outputFieldName = getFieldName(index);
    if (index == 0 && output.viskoresDataset.HasField(getFieldName(index, true))) {
        // if filter has created a dedicated output field, use it
        outputFieldName = getFieldName(index, true);
    }

    auto probe = std::make_unique<viskores::filter::resampling::Probe>();
    probe->SetGeometry(output.viskoresDataset);
    probe->SetOutputFieldName(outputFieldName);
    auto probeOutput = probe->Execute(input.viskoresDataset);

    // --------------------------------------------------------------------------------

    if (auto mapped = vtkmGetField(probeOutput, outputFieldName)) {
        std::cerr << "mapped data: " << *mapped << std::endl;
        updateMeta(mapped);
        mapped->copyAttributes(input.fields[index]);
        if (output.vistleGrid)
            mapped->setGrid(output.vistleGrid);
        output.fields[index] = mapped;
    } else {
        return Error("An error occurred while transforming the filter output field " + outputFieldName +
                     " to a Vistle object.");
    }

    return Success();
}

bool StreamlineVtkm::isValid(const ModuleStatusPtr &status) const
{
    if (strcmp(status->message(), ""))
        sendText(status->messageType(), status->message());

    return status->continueExecution();
}

bool StreamlineVtkm::tryToExecuteFilter(const std::unique_ptr<viskores::filter::Filter> &filter,
                                        const viskores::cont::DataSet &inputDataset,
                                        viskores::cont::DataSet &outputDataset) const
{
    std::string kind, description, message, backtrace;

    try {
        try {
            outputDataset = filter->Execute(inputDataset);
            return true;
        } catch (const viskores::cont::ErrorBadAllocation &error) {
            kind = "memory allocation error";
            description = "A memory allocation error occurred while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorBadDevice &error) {
            kind = "operation not supported by execution device";
            description = "The filter attempted to perform an operation that is not supported by the execution device";
            throw;
        } catch (const viskores::cont::ErrorBadType &error) {
            kind = "unsupported data type";
            description = "An unsupported data type was encountered while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorBadValue &error) {
            kind = "unsupported data value";
            description = "An invalid value was encountered while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorExecution &error) {
            kind = "execution environment error";
            description = "An error occurred in the execution environment while executing the filter";
            throw;
        } catch (const viskores::cont::ErrorFilterExecution &error) {
            kind = "filter setup error";
            description = "The filter has not been set up correctly";
            throw;
        } catch (const viskores::cont::ErrorInternal &error) {
            kind = "internal error";
            description = "An internal error occurred while executing the filter, indicating a bug in Viskores";
            throw;
        }
    } catch (const viskores::cont::Error &error) {
        if (kind.empty())
            kind = "unknown error";
        kind = "Viskores " + kind;
        message = error.GetMessage();
        backtrace = error.GetStackTrace();
    } catch (std::exception &error) {
        kind = "standard exception ";
        kind += typeid(error).name();
        message = error.what();
    } catch (...) {
        kind = "unknown exception";
    }

    std::stringstream msg;
    msg << "Execution of a Viskores filter failed with a " << kind << " exception";
    if (!description.empty())
        msg << " (" << description << ")";
    if (!message.empty())
        msg << ": " << message;
    sendError(msg.str());

    if (!backtrace.empty()) {
        msg << "\nBacktrace:\n" << backtrace;
    }

    std::cerr << msg.str() << std::endl;

    return false;
}
