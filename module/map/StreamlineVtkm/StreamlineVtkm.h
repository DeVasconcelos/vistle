#ifndef VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
#define VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H

#include <viskores/cont/ArrayHandle.h>
#include <viskores/Particle.h>

#include <vistle/vtkm/vtkm_module.h>

// TODO: find out why we get an error message when using more than 4x1x1 blocks
//       --> I assume because the Streamline filter can't execute on partitions in
//           concurrent threads seeing how CanThread() returns false...
// TODO: add max points?
// TODO: narrow conversion from vistle::Float to viskores::FloatDefault (or vistle::FLoat to Particles)
// TODO: calculate starting points on the GPU
// TODO: is it possible to also do backwards or bidirectional integration (like in Tracer)?
// TODO: I don't like that we have to convert input dataset + field again in prepareOutputField
class StreamlineVtkm: public vistle::VtkmModule {
public:
    StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm);

private:
    vistle::IntParameter *m_integrationMethod;
    vistle::IntParameter *m_numberOfPoints;
    vistle::IntParameter *m_numberOfSteps;
    vistle::IntParameter *m_startStyle;

    vistle::FloatParameter *m_stepSize;

    vistle::VectorParameter *m_direction;
    vistle::VectorParameter *m_startPoint1;
    vistle::VectorParameter *m_startPoint2;

    bool changeParameter(const vistle::Parameter *param) override;

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
