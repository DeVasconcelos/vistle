#ifndef VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
#define VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H

#include <vistle/vtkm/vtkm_module.h>

// TODO: find out why we get an error message when using more than 1x1x1 blocks
// TODO: add Line as StartStyle
// TODO: calculate starting points on the GPU
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

    ModuleStatusPtr prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                      const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                      viskores::cont::DataSet &dataset) const override;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;

    vistle::DataBase::ptr prepareOutputField(const viskores::cont::DataSet &dataset,
                                             const vistle::Object::const_ptr &inputGrid,
                                             const vistle::DataBase::const_ptr &inputField,
                                             const std::string &fieldName,
                                             const vistle::Object::const_ptr &outputGrid) const override;
};

#endif // VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
