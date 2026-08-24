#include <viskores/filter/flow/Streamline.h>
#include <viskores/filter/resampling/Probe.h>
#include <viskores/VectorAnalysis.h>

#include <vistle/util/enum.h>
#include <vistle/vtkm/convert.h>

#include "worklet/GenerateSeeds.h"

#include "StreamlineVtkm.h"

MODULE_MAIN(StreamlineVtkm)

using namespace vistle;

DEFINE_ENUM_WITH_STRING_CONVERSIONS(IntegrationMethod, (RK4)(Euler))
DEFINE_ENUM_WITH_STRING_CONVERSIONS(StartStyle, (Line)(Plane))

StreamlineVtkm::StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 3, MappedDataHandling::Require)
{
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

bool StreamlineVtkm::changeParameter(const Parameter *param)
{
    if (param == m_startStyle) {
        setParameterReadOnly(m_direction, m_startStyle->getValue() != Plane);
    } else if (param == m_maxNumberOfSeeds) {
        setParameterRange(m_numberOfSeeds, (Integer)1, m_maxNumberOfSeeds->getValue());
    }

    return Module::changeParameter(param);
}

ModuleStatusPtr StreamlineVtkm::prepareInputField(const Port *port, const Object::const_ptr &grid,
                                                  const DataBase::const_ptr &field, std::string &fieldName,
                                                  viskores::cont::DataSet &dataset) const
{
    if (port->getName() == "data_in") {
        if (auto in = Vec<Scalar, 3>::as(field))
            return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);

        return Error("Error: Input field must be a 3D vector field!");
    }

    return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);
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

DataBase::ptr StreamlineVtkm::prepareOutputField(const viskores::cont::DataSet &dataset,
                                                 const Object::const_ptr &inputGrid,
                                                 const DataBase::const_ptr &inputField, const std::string &fieldName,
                                                 const Object::const_ptr &outputGrid) const
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
