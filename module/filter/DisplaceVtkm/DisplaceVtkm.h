#ifndef VISTLE_DISPLACEVTKM_DISPLACEVTKM_H
#define VISTLE_DISPLACEVTKM_DISPLACEVTKM_H

#include <vistle/vtkm/vtkm_module.h>

class DisplaceVtkm: public vistle::VtkmModule {
public:
    DisplaceVtkm(const std::string &name, int moduleID, mpi::communicator comm);
    ~DisplaceVtkm();

private:
    vistle::IntParameter *p_component = nullptr;
    vistle::IntParameter *p_operation = nullptr;
    vistle::FloatParameter *p_scale = nullptr;

    ModuleStatusPtr prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                      const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                      viskores::cont::DataSet &dataset) const override;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;
};

#endif // VISTLE_DISPLACEVTKM_DISPLACEVTKM_H
