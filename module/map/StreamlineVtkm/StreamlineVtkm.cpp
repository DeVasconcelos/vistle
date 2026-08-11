#include <viskores/cont/ArrayHandle.h>
#include <viskores/filter/flow/Streamline.h>
#include <viskores/Particle.h>

#include "StreamlineVtkm.h"

MODULE_MAIN(StreamlineVtkm)

using namespace vistle;

StreamlineVtkm::StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 3, MappedDataHandling::Require)
{
    m_numberOfSteps = addIntParameter("number_of_steps", "number of integration steps", 100);
    m_stepSize = addFloatParameter("step_size", "integration step size", 0.1f);

    m_startPoint1 = addVectorParameter("startpoint1", "1st initial point", ParamVector(0, 0.2, 0));
    m_startPoint2 = addVectorParameter("startpoint2", "2nd initial point", ParamVector(1, 0, 0));
}

ModuleStatusPtr StreamlineVtkm::prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                                  const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                                  viskores::cont::DataSet &dataset) const
{
    if (auto in = Vec<Scalar, 3>::as(field))
        return VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);

    return Error("Error: Input field must be a 3D vector field!");
}

std::unique_ptr<viskores::filter::Filter> StreamlineVtkm::setUpFilter() const
{
    auto filter = std::make_unique<viskores::filter::flow::Streamline>();

    Vector3 startpoint1 = m_startPoint1->getValue();
    Vector3 startpoint2 = m_startPoint2->getValue();

    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(3);
    seedArray.WritePortal().Set(0, viskores::Particle({startpoint1[0], startpoint1[1], startpoint1[2]}, 0));
    seedArray.WritePortal().Set(1, viskores::Particle({startpoint2[0], startpoint2[1], startpoint2[2]}, 1));
    seedArray.WritePortal().Set(
        2, viskores::Particle({(startpoint2[0] - startpoint1[0]) / 2.f, (startpoint2[1] - startpoint1[1]) / 2.f,
                               (startpoint2[2] - startpoint1[2]) / 2.f},
                              2));

    filter->SetStepSize(m_stepSize->getValue());
    filter->SetNumberOfSteps(m_numberOfSteps->getValue());
    filter->SetSeeds(seedArray);

    return filter;
}
