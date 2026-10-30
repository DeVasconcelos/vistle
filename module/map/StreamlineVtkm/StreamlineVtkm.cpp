#include <algorithm>
#include <cmath>
#include <type_traits>

#include <viskores/TypeList.h>
#include <viskores/cont/Algorithm.h>
#include <viskores/cont/ArrayHandleStride.h>
#include <viskores/cont/EnvironmentTracker.h>
#include <viskores/cont/Invoker.h>
#include <viskores/cont/UnknownArrayHandle.h>
#include <viskores/filter/FieldSelection.h>
#include <viskores/filter/entity_extraction/Threshold.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/thirdparty/diy/diy.h>
#include <viskores/thirdparty/diy/mpi-cast.h>
#include <viskores/VectorAnalysis.h>
#include <viskores/worklet/WorkletMapField.h>

#include <vistle/util/enum.h>
#include <vistle/vtkm/convert.h>
#include <vistle/vtkm/vtkm_module_utils.h>

#include "worklet/GenerateSeeds.h"
#include "SampledStreamline.h"
#include "StreamlineVtkm.h"

using namespace vistle;

MODULE_MAIN(StreamlineVtkm)

DEFINE_ENUM_WITH_STRING_CONVERSIONS(IntegrationMethod, (Euler)(RK4))
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

bool StreamlineVtkm::GlobalData::getDataAtIndex(unsigned int index, std::vector<vistle::Object::const_ptr> &grids,
                                                std::vector<std::vector<vistle::DataBase::const_ptr>> &fieldsPerPort,
                                                viskores::cont::PartitionedDataSet &dataset)
{
    std::lock_guard<std::mutex> lock(this->mutex);
    if (index >= this->partitionedDatasets.size())
        return false;

    dataset = this->partitionedDatasets[index];
    if (dataset.GetNumberOfPartitions() == 0)
        return false;

    if (index < this->inputGrids.size())
        grids = this->inputGrids[index];

    fieldsPerPort.resize(this->inputFields.size());
    for (unsigned int port = 0; port < this->inputFields.size(); port++) {
        if (index < this->inputFields[port].size())
            fieldsPerPort[port] = this->inputFields[port][index];
    }

    return true;
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
: Module(name, moduleID, comm)
{
    createModulePorts();
    createModuleParameters();

    setReducePolicy(message::ReducePolicy::PerTimestep);
    setViskoresMPICommunicatorToCopyOf(comm);
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

bool StreamlineVtkm::compute(const std::shared_ptr<vistle::BlockTask> &task) const
{
    vistle::Object::const_ptr vistleGrid;
    std::vector<vistle::DataBase::const_ptr> fields;
    auto status = readInAndVerifyPorts(task, vistleGrid, fields);
    if (!checkAndNotify(status))
        return true;

    assert(m_outputPorts.size() == fields.size());

    auto timestep = vistleGrid->getTimestep();
    if (fields[0] && (timestep != fields[0]->getTimestep())) {
        sendError("timestep mismatch: grid = %d, field = %d", timestep, fields[0]->getTimestep());
        return true;
    }
    if (timestep < 0)
        timestep = -1;

    viskores::cont::DataSet viskoresDataset;
    status = transformInputToViskores(vistleGrid, fields, viskoresDataset);
    if (!checkAndNotify(status))
        return true;

    std::lock_guard<std::mutex> lock(m_globalData.mutex);
    auto numSteps = static_cast<std::size_t>(timestep + 1);
    m_globalData.resize(numSteps + 1, fields.size());

    m_globalData.partitionedDatasets[numSteps].AppendPartitions({viskoresDataset});
    m_globalData.inputGrids[numSteps].push_back(vistleGrid);
    for (std::size_t i = 0; i < fields.size(); ++i)
        m_globalData.inputFields[i][numSteps].push_back(fields[i]);

    return true;
}

bool StreamlineVtkm::reduce(int timestep)
{
    std::vector<vistle::Object::const_ptr> inputGrids;
    std::vector<std::vector<vistle::DataBase::const_ptr>> inputFields;
    viskores::cont::PartitionedDataSet inputPartitionedDataset;
    if (!m_globalData.getDataAtIndex(static_cast<unsigned int>(timestep + 1), inputGrids, inputFields,
                                     inputPartitionedDataset))
        return true;

    viskores::cont::PartitionedDataSet output;
    auto streamlineFilter = setUpFilter();
    if (!this->tryToExecuteFilter(*streamlineFilter, inputPartitionedDataset, output))
        return true;

    for (viskores::Id partitionIndex = 0; partitionIndex < output.GetNumberOfPartitions(); partitionIndex++) {
        auto dataset = output.GetPartition(partitionIndex);
        resolveMissingValues(dataset, inputPartitionedDataset);

        // streamlines with a single point (e.g., seeds that could not be advanced) would end up as spheres
        if (!removeDegenerateStreamlines(dataset))
            continue;
        auto outputGrid = vtkmGetGeometry(dataset);
        if (!outputGrid) {
            sendError("Could not convert StreamlineVtkm output geometry to a Vistle object.");
            continue;
        }

        outputGrid->setMeta(getMetaForOutput(timestep));
        updateMeta(outputGrid);
        // attributes (e.g. species name, color map range) are the same on every block, so let's just
        // use the first one (the Streamline filter does not tell us which input partition a streamline came from)
        if (!inputGrids.empty() && inputGrids.front())
            outputGrid->copyAttributes(inputGrids.front());

        // the filter samples the input fields on the input partition each streamline piece was traced in
        for (int port = 0; port < m_numPorts; port++) {
            if (!m_outputPorts[port]->isConnected())
                continue;

            DataBase::ptr field;
            if (dataset.HasPointField(getFieldName(port)))
                field = vtkmGetField(dataset, getFieldName(port));

            if (field) {
                field->setMeta(getMetaForOutput(timestep));
                updateMeta(field);
                // again, we simply take the first field to copy the attributes from, since we can't match the
                // input to the output partitions
                if (static_cast<std::size_t>(port) < inputFields.size() && !inputFields[port].empty()) {
                    if (auto inputField = inputFields[port].front()) {
                        auto mapping = field->mapping();
                        field->copyAttributes(inputField);
                        field->setMapping(mapping);
                    }
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

namespace {

// Calls f(T) with the floating point base component type of array and returns true, or returns false if array
// does not store floating point values (missing values are only marked for those, as NaN)
template<typename F>
bool withFloatType(const viskores::cont::UnknownArrayHandle &array, F &&f)
{
    if (array.IsBaseComponentType<viskores::Float32>()) {
        f(viskores::Float32());
        return true;
    }
    if (array.IsBaseComponentType<viskores::Float64>()) {
        f(viskores::Float64());
        return true;
    }
    return false;
}

struct MergeMissingWorklet: viskores::worklet::WorkletMapField {
    using ControlSignature = void(FieldIn srcHidden, FieldIn src, FieldInOut dst);
    using ExecutionSignature = void(_1, _2, _3);

    template<typename S, typename D>
    VISKORES_EXEC void operator()(viskores::UInt8 srcHidden, const S &src, D &dst) const
    {
        if (viskores::IsNan(dst) && !srcHidden)
            dst = static_cast<D>(src);
    }
};

// Copies the values of src to dst at all points where dst is NaN (i.e., could not be sampled while tracing) and
// where src is visible. Returns false if the arrays are incompatible.
bool mergeMissingValues(const viskores::cont::UnknownArrayHandle &dst, const viskores::cont::UnknownArrayHandle &src,
                        const viskores::cont::ArrayHandle<viskores::UInt8> &srcHidden)
{
    const viskores::Id n = dst.GetNumberOfValues();
    const viskores::IdComponent numComponents = dst.GetNumberOfComponentsFlat();
    if (src.GetNumberOfValues() != n || srcHidden.GetNumberOfValues() != n ||
        src.GetNumberOfComponentsFlat() != numComponents)
        return false;

    viskores::cont::Invoker invoke;
    bool called = false;
    withFloatType(dst, [&](auto d) {
        using D = decltype(d);
        auto tryType = [&](auto t) {
            using T = std::decay_t<decltype(t)>;
            if (called || !src.IsBaseComponentType<T>())
                return;
            called = true;
            for (viskores::IdComponent c = 0; c < numComponents; ++c) {
                // extracting without copying: the stride array aliases dst's buffer, so it is updated in place
                invoke(MergeMissingWorklet{}, srcHidden, src.ExtractComponent<T>(c),
                       dst.ExtractComponent<D>(c, viskores::CopyFlag::Off));
            }
        };
        viskores::ListForEach(tryType, viskores::TypeListScalarAll{});
    });
    return called;
}

struct IsMissingWorklet: viskores::worklet::WorkletMapField {
    using ControlSignature = void(FieldIn value, FieldOut missing);
    using ExecutionSignature = _2(_1);

    template<typename T>
    VISKORES_EXEC viskores::Id operator()(const T &value) const
    {
        return viskores::IsNan(value) ? 1 : 0;
    }
};

bool hasMissingValues(const viskores::cont::UnknownArrayHandle &array)
{
    viskores::cont::Invoker invoke;
    bool missing = false;
    withFloatType(array, [&](auto d) {
        using D = decltype(d);
        for (viskores::IdComponent c = 0; c < array.GetNumberOfComponentsFlat() && !missing; ++c) {
            viskores::cont::ArrayHandle<viskores::Id> isMissing;
            invoke(IsMissingWorklet{}, array.ExtractComponent<D>(c), isMissing);
            missing = viskores::cont::Algorithm::Reduce(isMissing, viskores::Id(0)) > 0;
            if (missing)
                break;
        }
    });
    return missing;
}

} // namespace

void StreamlineVtkm::resolveMissingValues(const viskores::cont::DataSet &streamlines,
                                          const viskores::cont::PartitionedDataSet &input) const
{
    // the filter samples each point on the block that traced it, which fails at block boundaries: probe those
    // points on all input blocks instead
    std::vector<std::string> names;
    for (auto i = 0; i < streamlines.GetNumberOfFields(); i++) {
        if (hasMissingValues(streamlines.GetField(i).GetData()))
            names.push_back(streamlines.GetField(i).GetName());
    }

    if (names.empty())
        return;

    // geometry only, so that the probed fields do not clash with the (incomplete) sampled ones
    viskores::cont::DataSet geometry;
    geometry.AddCoordinateSystem(streamlines.GetCoordinateSystem());
    geometry.SetCellSet(streamlines.GetCellSet());
    const auto geometryBounds = geometry.GetCoordinateSystem().GetBounds();

    for (viskores::Id block = 0; block < input.GetNumberOfPartitions(); ++block) {
        const auto &blockData = input.GetPartition(block);
        if (!blockData.GetCoordinateSystem().GetBounds().Intersection(geometryBounds).IsNonEmpty())
            continue;

        auto probeFilter = viskores::filter::resampling::Probe();
        probeFilter.SetGeometry(geometry);

        viskores::cont::DataSet probeOutput;
        if (!this->tryToExecuteFilter(probeFilter, blockData, probeOutput))
            continue;

        viskores::cont::ArrayHandle<viskores::UInt8> blockHidden;
        probeOutput.GetField("HIDDEN", viskores::cont::Field::Association::Points).GetData().AsArrayHandle(blockHidden);

        bool anyMissing = false;
        for (const auto &name: names) {
            if (!probeOutput.HasField(name, viskores::cont::Field::Association::Points))
                continue;
            auto dst = streamlines.GetPointField(name).GetData();
            mergeMissingValues(dst, probeOutput.GetField(name, viskores::cont::Field::Association::Points).GetData(),
                               blockHidden);
            anyMissing = anyMissing || hasMissingValues(dst);
        }
        if (!anyMissing)
            return;
    }
}

bool StreamlineVtkm::removeDegenerateStreamlines(viskores::cont::DataSet &streamlines) const
{
    const auto cellSet = streamlines.GetCellSet();
    const viskores::Id numCells = cellSet.GetNumberOfCells();

    viskores::cont::ArrayHandle<viskores::Float32> numPoints;
    numPoints.Allocate(numCells);
    auto portal = numPoints.WritePortal();
    bool anyDegenerate = false;
    for (viskores::Id cell = 0; cell < numCells; ++cell) {
        const auto n = cellSet.GetNumberOfPointsInCell(cell);
        portal.Set(cell, static_cast<viskores::Float32>(n));
        anyDegenerate = anyDegenerate || n < 2;
    }
    if (!anyDegenerate)
        return numCells > 0;

    const std::string name = "StreamlineNumPoints";
    streamlines.AddCellField(name, numPoints);
    viskores::filter::entity_extraction::Threshold threshold;
    threshold.SetActiveField(name, viskores::cont::Field::Association::Cells);
    threshold.SetThresholdAbove(2.0);
    viskores::cont::DataSet cleaned;
    if (!tryToExecuteFilter(threshold, streamlines, cleaned) || cleaned.GetNumberOfCells() == 0)
        return false;
    streamlines = cleaned;
    return true;
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

ModuleStatusPtr StreamlineVtkm::readInAndVerifyPorts(const std::shared_ptr<BlockTask> &task, Object::const_ptr &grid,
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

        // Streamlines filter expects data to be three-dimensional
        if (i == 0 && data && !Vec<Scalar, 3>::as(data))
            return Error("Input field at port " + m_inputPorts[i]->getName() + " must be a 3D vector field!");

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

ModuleStatusPtr StreamlineVtkm::transformInputToViskores(const Object::const_ptr &grid,
                                                         const std::vector<DataBase::const_ptr> &fields,
                                                         viskores::cont::DataSet &dataset) const
{
    auto status = vtkmSetGrid(dataset, grid);
    if (!status->continueExecution())
        return status;

    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (fields[i]) {
            status = vtkmAddField(dataset, fields[i], getFieldName(i));
            if (!status->continueExecution())
                return status;
        }
    }

    return Success();
}

namespace {
template<typename S>
viskores::Vec3f make_Vec3f(const vistle::ParameterVector<S> &v)
{
    typedef viskores::FloatDefault F;
    return viskores::Vec3f{static_cast<F>(v[0]), static_cast<F>(v[1]), static_cast<F>(v[2])};
}
} // namespace

viskores::cont::ArrayHandle<viskores::Particle> StreamlineVtkm::createSeedArray() const
{
    /*
        Since in Vistle VectorParameters (i.e., vector module parameters) are always vectors of
        doubles (independent of what the user sets VISTLE_DOUBLE_PRECISION to), and Viskores (in
        particular the array of viskores::Particles that the Streamline filter needs) strictly
        uses viskores::FloatDefault (which is dependent of VISTLE_DOUBLE_PRECISION), we have to
        cast the start points and direction to viskores::FloatDefault here.
    */
    viskores::Vec3f startpoint1 = make_Vec3f(m_startPoint1->getValue());
    viskores::Vec3f startpoint2 = make_Vec3f(m_startPoint2->getValue());

    if (m_startStyle->getValue() == StartStyle::Line) {
        return generateSeedsOnLine(m_numberOfSeeds->getValue(), startpoint1, startpoint2);

    } else {
        viskores::Vec3f direction = make_Vec3f(m_direction->getValue());
        return generateSeedsOnPlane(m_numberOfSeeds->getValue(), startpoint1, startpoint2, direction);
    }
}

std::unique_ptr<viskores::filter::Filter> StreamlineVtkm::setUpFilter() const
{
    auto filter = std::make_unique<vistle::SampledStreamline>();

    // sample the data of every connected output port along the streamlines
    std::vector<std::string> sampledFields;
    for (int port = 0; port < m_numPorts; ++port) {
        if (m_outputPorts[port]->isConnected())
            sampledFields.push_back(getFieldName(port));
    }
    filter->SetSampledFields(sampledFields);

    filter->SetStepSize(m_stepSize->getValue());
    filter->SetNumberOfSteps(m_numberOfSteps->getValue());

    auto seedArray = createSeedArray();
    filter->SetSeeds(seedArray);

    if (m_integrationMethod->getValue() == IntegrationMethod::Euler)
        filter->SetSolverEuler();
    else
        filter->SetSolverRK4();

    filter->SetActiveField(getFieldName(0));
    filter->SetOutputFieldName(getFieldName(0, true));
    filter->SetFieldsToPass("", viskores::cont::Field::Association::Any, viskores::filter::FieldSelection::Mode::All);

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

void StreamlineVtkm::createModulePorts()
{
    for (int i = 0; i < m_numPorts; ++i) {
        m_inputPorts.push_back(createInputPort(getPortName(i, false), "input grid with mapped data"));
        m_outputPorts.push_back(createOutputPort(getPortName(i, true), "output grid with mapped data"));

        linkPorts(m_inputPorts[i], m_outputPorts[i]);

        if (i > 0)
            setPortOptional(m_inputPorts[i], true);
    }
}

void StreamlineVtkm::createModuleParameters()
{
    /*
        Please note that the names of the module parameters need to match the parameter names
        of the Tracer module so that the two modules can be replaced with each other by right-
        clicking on them and that StreamlineVtkm can make use of COVER's TracerPlugin which
        allows the user to change the parameters in COVER.
    */
    setCurrentParameterGroup("Seed Points");
    const Integer max_no_seeds = 300;
    m_numberOfSeeds = addIntParameter("no_startp", "number of seed points", 2);
    setParameterRange(m_numberOfSeeds, (Integer)1, max_no_seeds);

    m_startStyle =
        addIntParameter("startStyle", "initial particle position configuration", StartStyle::Line, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(m_startStyle, StartStyle, );

    m_startPoint1 = addVectorParameter("startpoint1", "1st initial point", ParamVector(0, 0.2, 0));
    m_startPoint2 = addVectorParameter("startpoint2", "2nd initial point", ParamVector(1, 0, 0));

    m_direction = addVectorParameter("direction", "tracing direction", ParamVector(1, 0, 0));

    m_maxNumberOfSeeds =
        addIntParameter("max_no_startp", "maximum number of seeds (for parameter/slider limits)", max_no_seeds);
    setParameterRange(m_maxNumberOfSeeds, (Integer)2, (Integer)1000000);

    setCurrentParameterGroup("Stop Conditions");
    m_numberOfSteps = addIntParameter("steps_max", "maximum number of integration steps", 1000);

    setCurrentParameterGroup("Step Length Control");
    m_integrationMethod =
        addIntParameter("integration", "integration method", IntegrationMethod::RK4, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(m_integrationMethod, IntegrationMethod, );

    m_stepSize = addFloatParameter("h_init", "integration step size", 1e-03);
    setParameterRange("h_init", 0.0, 1e6);
}

void StreamlineVtkm::setViskoresMPICommunicatorToCopyOf(const mpi::communicator &comm) const
{
    MPI_Comm dupComm;
    MPI_Comm_dup(MPI_Comm(comm), &dupComm);
    viskoresdiy::mpi::communicator viskoresComm(viskoresdiy::mpi::make_DIY_MPI_Comm(dupComm), false);

    viskores::cont::EnvironmentTracker::SetCommunicator(viskoresComm);
}

std::string StreamlineVtkm::getFieldName(int index, bool output) const
{
    std::string name = "data_at_port_" + std::to_string(index);
    if (index == 0 && output)
        name += "_out";
    return name;
}

vistle::Meta StreamlineVtkm::getMetaForOutput(int timestep) const
{
    Meta meta;
    meta.setNumBlocks(size());
    meta.setBlock(rank());
    meta.setNumTimesteps(numTimesteps() > 0 ? numTimesteps() : -1);
    meta.setTimeStep(numTimesteps() > 0 ? timestep : -1);

    return meta;
}

std::string StreamlineVtkm::getPortName(int index, bool output) const
{
    std::string portName = output ? "data_out" : "data_in";
    portName += std::to_string(index);

    return portName;
}
