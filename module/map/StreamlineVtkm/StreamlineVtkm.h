#ifndef VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
#define VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H

#include <viskores/cont/ArrayHandle.h>
#include <viskores/Particle.h>

#include <vistle/vtkm/vtkm_module.h>

// TODO: narrow conversion from vistle::Float to viskores::FloatDefault (or vistle::Float to Particles)
// TODO: is it possible to also do backwards integration (like in Tracer)?
//       --> from what I can tell, the filter only supports forward integration as neither the filter itself
//           nor the integrator worklets nor the Particle data structure itself have an option for backwards
//           integration. So we would have to implement that ourselves.
class StreamlineVtkm: public vistle::VtkmModule {
public:
    StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm);

private:
    vistle::IntParameter *m_integrationMethod;
    vistle::IntParameter *m_numberOfSeeds, *m_maxNumberOfSeeds, *m_numberOfSteps;
    vistle::IntParameter *m_startStyle;

    vistle::FloatParameter *m_stepSize;

    vistle::VectorParameter *m_direction;
    vistle::VectorParameter *m_startPoint1, *m_startPoint2;

    bool changeParameter(const vistle::Parameter *param) override;

    ModuleStatusPtr prepareInputGrid(const vistle::Object::const_ptr &grid,
                                     viskores::cont::DataSet &dataset) const override;

    ModuleStatusPtr prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                      const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                      viskores::cont::DataSet &dataset) const override;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;

    vistle::DataBase::ptr prepareOutputField(const viskores::cont::DataSet &dataset,
                                             const vistle::Object::const_ptr &inputGrid,
                                             const vistle::DataBase::const_ptr &inputField,
                                             const std::string &fieldName,
                                             const vistle::Object::const_ptr &outputGrid) const override;

    viskores::cont::ArrayHandle<viskores::Particle> createSeedArray() const;
};

#endif // VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
