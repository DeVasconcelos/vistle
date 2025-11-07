#ifndef VISTLE_DISPLACEVTKM_DISPLACEVTKM_H
#define VISTLE_DISPLACEVTKM_DISPLACEVTKM_H

#include <vistle/vtkm/vtkm_module.h>

// TODO: make this a regular module (or better yet, just add it to MapDrape)
// we need to work on a non-const dataset, but can't do it with the VtkmModule design
// (the code in prepareInputField can be reused)
class DisplaceVtkm: public vistle::VtkmModule {
public:
    DEFINE_ENUM_WITH_STRING_CONVERSIONS(DisplaceComponent, (X)(Y)(Z)(All))
    DEFINE_ENUM_WITH_STRING_CONVERSIONS(DisplaceOperation, (Set)(Add)(Multiply))

    DisplaceVtkm(const std::string &name, int moduleID, mpi::communicator comm);
    ~DisplaceVtkm();

private:
    vistle::IntParameter *p_component = nullptr;
    vistle::IntParameter *p_operation = nullptr;
    vistle::FloatParameter *p_scale = nullptr;

    template<typename ScalarArrayType, typename CoordsArrayType>
    void applyOperation(const ScalarArrayType &scalar, CoordsArrayType &coords, viskores::IdComponent c,
                        DisplaceOperation operation, viskores::FloatDefault scale) const;

    ModuleStatusPtr prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                      const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                      viskores::cont::DataSet &dataset) const override;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;
};

#endif // VISTLE_DISPLACEVTKM_DISPLACEVTKM_H
