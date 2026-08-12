#include <viskores/cont/ArrayHandle.h>
#include <viskores/filter/flow/Streamline.h>
#include <viskores/Particle.h>

#include "StreamlineVtkm.h"

MODULE_MAIN(StreamlineVtkm)

using namespace vistle;

StreamlineVtkm::StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 3, MappedDataHandling::Require)
{
    m_numberOfPoints = addIntParameter("number_of_points", "number of seed points", 4);
    m_numberOfSteps = addIntParameter("number_of_steps", "number of integration steps", 100);
    m_stepSize = addFloatParameter("step_size", "integration step size", 0.1f);

    m_direction = addVectorParameter("direction", "tracing direction", ParamVector(0, 0, 1));
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

std::vector<Vector3> calculateStartingPoint(Index numpoints, const Vector3 &startpoint1, const Vector3 &startpoint2,
                                            const Vector3 &direction)
{
    std::vector<Vector3> startpoints;
    auto normedDirection = direction;
    normedDirection.normalize();
    Vector3 v = startpoint2 - startpoint1;
    Scalar l = normedDirection.dot(v);
    Vector3 v0 = normedDirection * l;
    Vector3 v1 = v - v0;
    Scalar r = v0.norm();
    Scalar s = v1.norm();
    Index n0, n1;
    if (r > s) {
        n1 = Index(sqrt(numpoints * s / r)) + 1;
        if (n1 <= 1)
            n1 = 2;
        n0 = numpoints / n1;
        if (n0 <= 1)
            n0 = 2;
    } else {
        n0 = Index(sqrt(numpoints * r / s)) + 1;
        if (n0 <= 1)
            n0 = 2;
        n1 = numpoints / n0;
        if (n1 <= 1)
            n1 = 2;
    }

    numpoints = n0 * n1;
    startpoints.resize(numpoints);

    Scalar s0 = Scalar(1) / (n0 - 1);
    Scalar s1 = Scalar(1) / (n1 - 1);
    for (Index i = 0; i < n0; ++i) {
        for (Index j = 0; j < n1; ++j) {
            startpoints[i * n1 + j] = startpoint1 + v0 * s0 * i + v1 * s1 * j;
        }
    }

    return startpoints;
}

std::unique_ptr<viskores::filter::Filter> StreamlineVtkm::setUpFilter() const
{
    auto filter = std::make_unique<viskores::filter::flow::Streamline>();

    auto numPoints = m_numberOfPoints->getValue();
    auto points = calculateStartingPoint(numPoints, m_startPoint1->getValue(), m_startPoint2->getValue(),
                                         m_direction->getValue());


    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(numPoints);
    for (Index i = 0; i < numPoints; i++)
        seedArray.WritePortal().Set(i, viskores::Particle({points[i][0], points[i][1], points[i][2]}, i));

    filter->SetStepSize(m_stepSize->getValue());
    filter->SetNumberOfSteps(m_numberOfSteps->getValue());
    filter->SetSeeds(seedArray);

    return filter;
}

vistle::DataBase::ptr StreamlineVtkm::prepareOutputField(const viskores::cont::DataSet &dataset,
                                                         const vistle::Object::const_ptr &inputGrid,
                                                         const vistle::DataBase::const_ptr &inputField,
                                                         const std::string &fieldName,
                                                         const vistle::Object::const_ptr &outputGrid) const
{
    return nullptr;
}
