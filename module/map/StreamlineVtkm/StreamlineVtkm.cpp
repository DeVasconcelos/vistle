#include <viskores/filter/flow/Streamline.h>
#include <viskores/filter/resampling/Probe.h>

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

viskores::cont::ArrayHandle<viskores::Particle> generateSeedsOnLine(Index numSeeds, const Vector3 &startpoint1,
                                                                    const Vector3 &startpoint2)
{
    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(numSeeds);

    if (numSeeds == 1) {
        auto point = (startpoint1 + startpoint2) * 0.5;
        seedArray.WritePortal().Set(0, viskores::Particle({point[0], point[1], point[2]}, 0));
    } else {
        Vector3 delta = (startpoint2 - startpoint1) / (numSeeds - 1);
        for (Index i = 0; i < numSeeds; i++) {
            auto point = startpoint1 + i * delta;
            seedArray.WritePortal().Set(i, viskores::Particle({point[0], point[1], point[2]}, i));
        }
    }

    return seedArray;
}

void computePlaneSpans(const Vector3 &startpoint1, const Vector3 &startpoint2, const Vector3 &direction,
                       Vector3 &parallelSpan, Vector3 &orthogonalSpan)
{
    auto normedDirection = direction;
    normedDirection.normalize();
    Vector3 seedSpan = startpoint2 - startpoint1;
    parallelSpan = normedDirection * normedDirection.dot(seedSpan);
    orthogonalSpan = seedSpan - parallelSpan;
}

void computeSeedCountPerSpan(Index numSeeds, Scalar length0, Scalar length1, Index &count0, Index &count1)
{
    if (length0 > length1) {
        count1 = Index(sqrt(numSeeds * length1 / length0)) + 1;
        if (count1 <= 1)
            count1 = 2;
        count0 = numSeeds / count1;
        if (count0 <= 1)
            count0 = 2;
    } else {
        count0 = Index(sqrt(numSeeds * length0 / length1)) + 1;
        if (count0 <= 1)
            count0 = 2;
        count1 = numSeeds / count0;
        if (count1 <= 1)
            count1 = 2;
    }
}

viskores::cont::ArrayHandle<viskores::Particle>
generateSeedsOnPlane(Index numSeeds, const Vector3 &startpoint1, const Vector3 &startpoint2, const Vector3 &direction)
{
    Vector3 parallelSpan, orthogonalSpan;
    computePlaneSpans(startpoint1, startpoint2, direction, parallelSpan, orthogonalSpan);

    Index n0, n1;
    computeSeedCountPerSpan(numSeeds, parallelSpan.norm(), orthogonalSpan.norm(), n0, n1);

    numSeeds = n0 * n1;
    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(numSeeds);

    Scalar s0 = Scalar(1) / (n0 - 1);
    Scalar s1 = Scalar(1) / (n1 - 1);
    for (Index i = 0; i < n0; ++i) {
        for (Index j = 0; j < n1; ++j) {
            auto point = startpoint1 + parallelSpan * s0 * i + orthogonalSpan * s1 * j;
            seedArray.WritePortal().Set(i * n1 + j, viskores::Particle({point[0], point[1], point[2]}, i * n1 + j));
        }
    }

    return seedArray;
}

viskores::cont::ArrayHandle<viskores::Particle> StreamlineVtkm::createSeedArray() const
{
    if (m_startStyle->getValue() == StartStyle::Line)
        return generateSeedsOnLine(m_numberOfPoints->getValue(), m_startPoint1->getValue(), m_startPoint2->getValue());
    else
        return generateSeedsOnPlane(m_numberOfPoints->getValue(), m_startPoint1->getValue(), m_startPoint2->getValue(),
                                    m_direction->getValue());
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
