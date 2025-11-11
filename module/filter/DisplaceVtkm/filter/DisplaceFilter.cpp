#include <string>
#include <type_traits>

#include <viskores/cont/ArrayHandleExtractComponent.h>
#include <viskores/cont/ArrayHandleSOA.h>

#include "DisplaceFilter.h"
#include "DisplaceWorklet.h"

VISKORES_CONT DisplaceFilter::DisplaceFilter()
: m_component(DisplaceComponent::X), m_operation(DisplaceOperation::Add), m_scale(1.0f)
{}

// applies desired displace operation to `coords`. It is assumed that `field` and `coords` have the same dimension.
template<typename FieldArrayType, typename CoordsArrayType>
void applyDisplaceOperation(const FieldArrayType &field, CoordsArrayType &coords, viskores::FloatDefault scale,
                            DisplaceFilter::DisplaceOperation operation)
{
    static_assert(!std::is_const<std::remove_reference_t<CoordsArrayType>>::value,
                  "CoordsArrayType must be non-const as its contents will be modified by applyDisplaceOperation!");

    viskores::cont::Invoker invoke;
    switch (operation) {
    case DisplaceFilter::DisplaceOperation::Set:
        invoke(SetDisplaceWorklet{scale}, field, coords, coords);
        break;
    case DisplaceFilter::DisplaceOperation::Add:
        invoke(AddDisplaceWorklet{scale}, field, coords, coords);
        break;
    case DisplaceFilter::DisplaceOperation::Multiply:
        invoke(MultiplyDisplaceWorklet{scale}, field, coords, coords);
        break;
    default:
        throw viskores::cont::ErrorBadValue("Error in DisplaceVtkm: Encountered unknown DisplaceOperation value!");
    }
}

// applies desired displace operation to a specific component of the `coords` array
template<typename FieldArrayType, typename CoordsArrayType>
void applyDisplaceOperation(const FieldArrayType &field, CoordsArrayType &coords, viskores::FloatDefault scale,
                            DisplaceFilter::DisplaceOperation operation, viskores::IdComponent c)
{
    auto desiredComponent = viskores::cont::ArrayHandleExtractComponent<CoordsArrayType>(coords, c);

    applyDisplaceOperation(field, desiredComponent, scale, operation);
}

VISKORES_CONT viskores::cont::DataSet DisplaceFilter::DoExecute(const viskores::cont::DataSet &inputDataset)
{
    auto inputField = this->GetFieldFromDataSet(inputDataset);
    auto inputCoords = inputDataset.GetCoordinateSystem(this->GetActiveCoordinateSystemIndex());

    viskores::cont::UnknownArrayHandle outputCoords;

    // we assume the point coordinates are three dimensional
    constexpr viskores::IdComponent COORDS_DIM = 3;
    this->CastAndCallVecField<COORDS_DIM>(inputCoords.GetData(), [&](const auto &coords) {
        using CoordsArrayType = std::decay_t<decltype(coords)>;
        using CoordType = typename CoordsArrayType::ValueType;

        // we assume the data field is either 1D or 3D
        using TypeListMappedData =
            viskores::List<viskores::Float32, viskores::Float64, viskores::Vec3f_32, viskores::Vec3f_64>;
        inputField.GetData().CastAndCallForTypesWithFloatFallback<TypeListMappedData, VISKORES_DEFAULT_STORAGE_LIST>(
            [&](const auto &field) {
                using FieldType = typename std::decay_t<decltype(field)>::ValueType;
                constexpr bool fieldIsScalar = std::is_arithmetic_v<FieldType>;

                viskores::cont::ArrayHandle<CoordType> result;
                result.Allocate(coords.GetNumberOfValues());
                viskores::cont::ArrayCopy(coords, result);

                if constexpr (fieldIsScalar) {
                    if (m_component == DisplaceComponent::All) {
                        for (viskores::IdComponent c = 0; c < COORDS_DIM; c++)
                            applyDisplaceOperation(field, result, m_scale, m_operation, c);

                        outputCoords = result;
                    } else if (m_component == DisplaceComponent::X || m_component == DisplaceComponent::Y ||
                               m_component == DisplaceComponent::Z) {
                        viskores::IdComponent c = static_cast<viskores::IdComponent>(m_component);
                        applyDisplaceOperation(field, result, m_scale, m_operation, c);

                        outputCoords = result;

                    } else {
                        throw viskores::cont::ErrorBadValue(
                            "Error in DisplaceFilter: Encountered unknown DisplaceComponent value!");
                    }
                } else if constexpr (!fieldIsScalar && FieldType::NUM_COMPONENTS == COORDS_DIM) {
                    applyDisplaceOperation(field, result, m_scale, m_operation);

                    outputCoords = result;
                } else {
                    throw viskores::cont::ErrorBadValue(
                        "Error in DisplaceFilter: Cannot apply filter on point coordinates of dimension " +
                        std::to_string(COORDS_DIM) + " with data field of dimension " +
                        std::to_string(FieldType::NUM_COMPONENTS) + "!");
                }
            });
    });

    auto outputFieldName = this->GetOutputFieldName();
    if (outputFieldName == "")
        outputFieldName = inputField.GetName() + "_displaced";

    return this->CreateResultCoordinateSystem(
        inputDataset, inputDataset.GetCellSet(), inputCoords.GetName(), outputCoords,
        [&](viskores::cont::DataSet &out, const viskores::cont::Field &fieldToPass) {
            out.AddField(viskores::cont::Field(
                fieldToPass.GetName() == this->GetActiveFieldName() ? outputFieldName : fieldToPass.GetName(),
                fieldToPass.GetAssociation(), fieldToPass.GetData()));
        });
}
