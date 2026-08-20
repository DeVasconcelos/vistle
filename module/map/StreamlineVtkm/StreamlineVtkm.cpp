#include <viskores/filter/flow/Streamline.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/VectorAnalysis.h>
#include <viskores/worklet/DispatcherMapField.h>
#include <viskores/worklet/WorkletMapField.h>

#include <vistle/util/enum.h>
#include <vistle/vtkm/convert.h>

#include "StreamlineVtkm.h"

MODULE_MAIN(StreamlineVtkm)

using namespace vistle;

DEFINE_ENUM_WITH_STRING_CONVERSIONS(IntegrationMethod, (RK4)(Euler))
DEFINE_ENUM_WITH_STRING_CONVERSIONS(StartStyle, (Line)(Plane))

StreamlineVtkm::StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 3, MappedDataHandling::Require)
{
    setCurrentParameterGroup("Seed Points");
    m_numberOfPoints = addIntParameter("number_of_seeds", "number of seed points", 2);
    m_startStyle =
        addIntParameter("start_style", "initial particle position configuration", StartStyle::Line, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(m_startStyle, StartStyle, );
    m_startPoint1 = addVectorParameter("startpoint1", "1st initial point", ParamVector(0, 0.2, 0));
    m_startPoint2 = addVectorParameter("startpoint2", "2nd initial point", ParamVector(1, 0, 0));
    m_direction = addVectorParameter("direction", "tracing direction", ParamVector(0, 0, 1));

    setCurrentParameterGroup("Stop Conditions");
    m_numberOfSteps = addIntParameter("steps_max", "maximum number of integration steps", 100);

    setCurrentParameterGroup("Step Length Control");
    m_integrationMethod =
        addIntParameter("integration_method", "integration method", IntegrationMethod::RK4, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(m_integrationMethod, IntegrationMethod, );
    m_stepSize = addFloatParameter("step_size", "integration step size", 0.1f);
}

bool StreamlineVtkm::changeParameter(const vistle::Parameter *param)
{
    if (param == m_startStyle)
        setParameterReadOnly(m_direction, m_startStyle->getValue() != Plane);

    return Module::changeParameter(param);
}

ModuleStatusPtr StreamlineVtkm::prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                                  const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                                  viskores::cont::DataSet &dataset) const
{
    if (port->getName() == "data_in") {
        if (auto in = Vec<Scalar, 3>::as(field))
            return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);

        return Error("Error: Input field must be a 3D vector field!");
    }

    return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);
}

class GenerateSeedsOnLineWorklet: public viskores::worklet::WorkletMapField {
public:
    VISKORES_CONT GenerateSeedsOnLineWorklet(viskores::Id numSeeds, const viskores::Vec3f &startpoint1,
                                             const viskores::Vec3f &startpoint2)
    : m_numSeeds(numSeeds)
    , m_startpoint1(startpoint1)
    , m_startpoint2(startpoint2)
    , m_delta((m_startpoint2 - m_startpoint1) / (m_numSeeds - 1))
    {}

    VISKORES_CONT GenerateSeedsOnLineWorklet(viskores::Id numSeeds, const viskores::Vec3f &startpoint1,
                                             const viskores::Vec3f &startpoint2, const viskores::Vec3f &delta)
    : m_numSeeds(numSeeds), m_startpoint1(startpoint1), m_startpoint2(startpoint2), m_delta(delta)
    {}

    using ControlSignature = void(FieldOut seeds);
    using ExecutionSignature = void(_1, WorkIndex);

    VISKORES_EXEC void operator()(viskores::Particle &seed, const viskores::Id &index) const
    {
        seed = viskores::Particle({m_startpoint1[0] + index * m_delta[0], m_startpoint1[1] + index * m_delta[1],
                                   m_startpoint1[2] + index * m_delta[2]},
                                  index);
    }

private:
    viskores::Id m_numSeeds;

    viskores::Vec3f m_startpoint1;
    viskores::Vec3f m_startpoint2;
    viskores::Vec3f m_delta;
};

viskores::cont::ArrayHandle<viskores::Particle>
generateSeedsOnLine(viskores::Id numSeeds, const viskores::Vec3f &startpoint1, const viskores::Vec3f &startpoint2)
{
    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(numSeeds);

    if (numSeeds == 1) {
        auto point = (startpoint1 + startpoint2) * 0.5;
        seedArray.WritePortal().Set(0, viskores::Particle({point[0], point[1], point[2]}, 0));
    } else {
        viskores::worklet::DispatcherMapField<GenerateSeedsOnLineWorklet>(
            GenerateSeedsOnLineWorklet(numSeeds, startpoint1, startpoint2))
            .Invoke(seedArray);
    }

    return seedArray;
}

void computePlaneSpans(const viskores::Vec3f &startpoint1, const viskores::Vec3f &startpoint2,
                       const viskores::Vec3f &direction, viskores::Vec3f &parallelSpan, viskores::Vec3f &orthogonalSpan)
{
    auto normedDirection = direction;
    viskores::Normalize(normedDirection);
    auto seedSpan = startpoint2 - startpoint1;
    parallelSpan = normedDirection * viskores::Dot(normedDirection, seedSpan);
    orthogonalSpan = seedSpan - parallelSpan;
}

void computeSeedCountPerSpan(viskores::Id numSeeds, viskores::FloatDefault length0, viskores::FloatDefault length1,
                             viskores::Id &count0, viskores::Id &count1)
{
    if (length0 > length1) {
        count1 = viskores::Id(viskores::Sqrt(numSeeds * length1 / length0)) + 1;
        if (count1 <= 1)
            count1 = 2;
        count0 = numSeeds / count1;
        if (count0 <= 1)
            count0 = 2;
    } else {
        count0 = viskores::Id(viskores::Sqrt(numSeeds * length0 / length1)) + 1;
        if (count0 <= 1)
            count0 = 2;
        count1 = numSeeds / count0;
        if (count1 <= 1)
            count1 = 2;
    }
}

class GenerateSeedsOnPlaneWorklet: public viskores::worklet::WorkletMapField {
public:
    VISKORES_CONT
    GenerateSeedsOnPlaneWorklet(viskores::Id numSeeds, viskores::Id n0, viskores::Id n1,
                                const viskores::Vec3f &startpoint1, const viskores::Vec3f &startpoint2,
                                const viskores::Vec3f &parallelSpan, const viskores::Vec3f &orthogonalSpan)
    : m_numSeeds(numSeeds)
    , m_n0(n0)
    , m_n1(n1)
    , m_startpoint1(startpoint1)
    , m_startpoint2(startpoint2)
    , m_parallelSpan(parallelSpan)
    , m_orthogonalSpan(orthogonalSpan)
    {}

    using ControlSignature = void(FieldOut seeds);
    using ExecutionSignature = void(_1, WorkIndex);

    VISKORES_EXEC void operator()(viskores::Particle &seed, const viskores::Id &index) const
    {
        viskores::Id i = index / m_n1;
        viskores::Id j = index % m_n1;

        Scalar s0 = Scalar(1) / (m_n0 - 1);
        Scalar s1 = Scalar(1) / (m_n1 - 1);
        auto point = m_startpoint1 + m_parallelSpan * s0 * i + m_orthogonalSpan * s1 * j;
        seed = viskores::Particle({point[0], point[1], point[2]}, index);
    }

private:
    viskores::Id m_numSeeds;

    viskores::Id m_n0;
    viskores::Id m_n1;

    viskores::Vec3f m_startpoint1;
    viskores::Vec3f m_startpoint2;
    viskores::Vec3f m_parallelSpan;
    viskores::Vec3f m_orthogonalSpan;
};

viskores::cont::ArrayHandle<viskores::Particle> generateSeedsOnPlane(Index numSeeds, const viskores::Vec3f &startpoint1,
                                                                     const viskores::Vec3f &startpoint2,
                                                                     const viskores::Vec3f &direction)
{
    viskores::Vec3f parallelSpan, orthogonalSpan;
    computePlaneSpans(startpoint1, startpoint2, direction, parallelSpan, orthogonalSpan);

    viskores::Id n0, n1;
    computeSeedCountPerSpan(numSeeds, viskores::Magnitude(parallelSpan), viskores::Magnitude(orthogonalSpan), n0, n1);

    numSeeds = n0 * n1;
    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(numSeeds);

    viskores::worklet::DispatcherMapField<GenerateSeedsOnPlaneWorklet>(
        GenerateSeedsOnPlaneWorklet(numSeeds, n0, n1, startpoint1, startpoint2, parallelSpan, orthogonalSpan))
        .Invoke(seedArray);

    return seedArray;
}

viskores::cont::ArrayHandle<viskores::Particle> StreamlineVtkm::createSeedArray() const
{
    viskores::Vec3f startpoint1{m_startPoint1->getValue()[0], m_startPoint1->getValue()[1],
                                m_startPoint1->getValue()[2]};
    viskores::Vec3f startpoint2{m_startPoint2->getValue()[0], m_startPoint2->getValue()[1],
                                m_startPoint2->getValue()[2]};

    if (m_startStyle->getValue() == StartStyle::Line) {
        return generateSeedsOnLine(m_numberOfPoints->getValue(), startpoint1, startpoint2);

    } else {
        viskores::Vec3f direction{m_direction->getValue()[0], m_direction->getValue()[1], m_direction->getValue()[2]};
        return generateSeedsOnPlane(m_numberOfPoints->getValue(), startpoint1, startpoint2, direction);
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

vistle::DataBase::ptr StreamlineVtkm::prepareOutputField(const viskores::cont::DataSet &dataset,
                                                         const vistle::Object::const_ptr &inputGrid,
                                                         const vistle::DataBase::const_ptr &inputField,
                                                         const std::string &fieldName,
                                                         const vistle::Object::const_ptr &outputGrid) const
{
    // The Streamline filter only returns a geometry of polylines. To match the Tracer behavior, we
    // need to resample the input field onto the output geometry. To keep this on the device, we use
    // the Probe filter for this.
    // FIXME: Due to the rigid structure of VtkmModule, we need to re-convert the input grid and field
    // to a Viskores dataset again here (this already happens in the prepareInput-methods).
    viskores::cont::DataSet inputDataset;
    auto status = vtkmSetGrid(inputDataset, inputGrid);
    if (!isValid(status))
        return nullptr;

    status = vtkmAddField(inputDataset, inputField, fieldName);
    if (!isValid(status))
        return nullptr;

    auto probe = std::make_unique<viskores::filter::resampling::Probe>();
    probe->SetGeometry(dataset);
    probe->SetOutputFieldName(fieldName);
    auto probeOutput = probe->Execute(inputDataset);

    return VtkmModule::prepareOutputField(probeOutput, inputGrid, inputField, fieldName, outputGrid);
}
