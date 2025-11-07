#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/ArrayHandleExtractComponent.h>
#include <viskores/cont/Invoker.h>

#include "filter/DisplaceFilter.h"
#include "filter/DisplaceWorklet.h"

#include "DisplaceVtkm.h"

MODULE_MAIN(DisplaceVtkm)

using namespace vistle;

DisplaceVtkm::DisplaceVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 5, MappedDataHandling::Require)
{
    p_component = addIntParameter("component", "component to displace for scalar input",
                                  DisplaceFilter::DisplaceComponent::Z, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(p_component, DisplaceComponent, DisplaceFilter);

    p_operation = addIntParameter("operation", "displacement operation to apply to selected component or element-wise",
                                  DisplaceFilter::DisplaceOperation::Add, Parameter::Choice);
    V_ENUM_SET_CHOICES_SCOPE(p_operation, DisplaceOperation, DisplaceFilter);

    p_scale = addFloatParameter("scale", "scaling factor for displacement", 1.);
}

DisplaceVtkm::~DisplaceVtkm()
{}

template<viskores::IdComponent VecSize>
struct ScalarToVec {
    template<typename T>
    using type = viskores::Vec<T, VecSize>;
};

template<typename ScalarArrayType, typename CoordsArrayType>
void applyOperation(const ScalarArrayType &scalar, CoordsArrayType &coords, viskores::FloatDefault scale,
                    DisplaceFilter::DisplaceOperation operation, viskores::IdComponent c)
{
    auto desiredComponent = viskores::cont::ArrayHandleExtractComponent<CoordsArrayType>(coords, c);
    using ComponentType = typename std::decay_t<decltype(desiredComponent)>::ValueType;


    viskores::cont::ArrayHandle<ComponentType> result;

    viskores::cont::Invoker invoke;
    switch (operation) {
    case DisplaceFilter::DisplaceOperation::Set:
        invoke(SetDisplaceWorklet{scale}, scalar, desiredComponent, result);
        break;
    case DisplaceFilter::DisplaceOperation::Add:
        invoke(AddDisplaceWorklet{scale}, scalar, desiredComponent, result);
        break;
    case DisplaceFilter::DisplaceOperation::Multiply:
        invoke(MultiplyDisplaceWorklet{scale}, scalar, desiredComponent, result);
        break;
    default:
        throw viskores::cont::ErrorBadValue("Error in DisplaceVtkm: Encountered unknown DisplaceOperation value!");
    }

    viskores::cont::ArrayCopy(result, desiredComponent);
}

ModuleStatusPtr DisplaceVtkm::prepareInputField(const vistle::Port *port, const vistle::Object::const_ptr &grid,
                                                const vistle::DataBase::const_ptr &field, std::string &fieldName,
                                                viskores::cont::DataSet &dataset) const
{
    auto status = VtkmModule::prepareInputField(port, grid, field, fieldName, dataset);
    if (!isValid(status))
        return status;


    if (fieldName == getFieldName(0)) {
        auto component = static_cast<DisplaceFilter::DisplaceComponent>(p_component->getValue());
        auto operation = static_cast<DisplaceFilter::DisplaceOperation>(p_operation->getValue());
        auto scale = static_cast<viskores::FloatDefault>(p_scale->getValue());


        auto scalarField = dataset.GetField(fieldName).GetData();
        auto coords = dataset.GetCoordinateSystem().GetData();

        try {
            scalarField.CastAndCallForTypesWithFloatFallback<viskores::TypeListFieldScalar,
                                                             VISKORES_DEFAULT_STORAGE_LIST>([&](const auto &scalar) {
                // TODO: don't hard code...
                constexpr viskores::Id N = 3;
                using VecList = viskores::ListTransform<viskores::TypeListFieldScalar, ScalarToVec<N>::template type>;
                coords.CastAndCallForTypesWithFloatFallback<VecList, VISKORES_DEFAULT_STORAGE_LIST>([&](auto &coords) {
                    if (component == DisplaceFilter::DisplaceComponent::All) {
                        for (viskores::IdComponent c = 0; c < N; c++)
                            applyOperation(scalar, coords, scale, operation, c);

                    } else if (component == DisplaceFilter::DisplaceComponent::X ||
                               component == DisplaceFilter::DisplaceComponent::Y ||
                               component == DisplaceFilter::DisplaceComponent::Z) {
                        applyOperation(scalar, coords, scale, operation,
                                       static_cast<viskores::IdComponent>(p_component->getValue()));

                    } else {
                        throw viskores::cont::ErrorBadValue(
                            "Error in DisplaceVtkm: Encountered unknown DisplaceComponent value!");
                    }
                });
            });
        } catch (const viskores::cont::Error &error) {
            return Error(error.GetMessage());
        }
    }
    return Success();
}

std::unique_ptr<viskores::filter::Filter> DisplaceVtkm::setUpFilter() const
{
    return nullptr;
}
