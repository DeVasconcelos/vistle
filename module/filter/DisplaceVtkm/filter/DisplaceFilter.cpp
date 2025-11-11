#include <string>
#include <type_traits>

#include <viskores/cont/ArrayHandleExtractComponent.h>
#include <viskores/cont/ArrayHandleSOA.h>

#include "DisplaceFilter.h"
#include "DisplaceWorklet.h"

VISKORES_CONT DisplaceFilter::DisplaceFilter()
: m_component(DisplaceComponent::X), m_operation(DisplaceOperation::Add), m_scale(1.0f)
{}

template<typename ScalarArrayType, typename CoordsArrayType>
void applyOperation(const ScalarArrayType &scalar, CoordsArrayType &coords, viskores::FloatDefault scale,
                    DisplaceFilter::DisplaceOperation operation, viskores::IdComponent c)
{
    static_assert(!std::is_const<std::remove_reference_t<CoordsArrayType>>::value,
                  "CoordsArrayType must be non-const as its contents will be modified by applyOperation!");

    auto desiredComponent = viskores::cont::ArrayHandleExtractComponent<CoordsArrayType>(coords, c);

    viskores::cont::Invoker invoke;
    switch (operation) {
    case DisplaceFilter::DisplaceOperation::Set:
        invoke(SetDisplaceWorklet{scale}, scalar, desiredComponent, desiredComponent);
        break;
    case DisplaceFilter::DisplaceOperation::Add:
        invoke(AddDisplaceWorklet{scale}, scalar, desiredComponent, desiredComponent);
        break;
    case DisplaceFilter::DisplaceOperation::Multiply:
        invoke(MultiplyDisplaceWorklet{scale}, scalar, desiredComponent, desiredComponent);
        break;
    default:
        throw viskores::cont::ErrorBadValue("Error in DisplaceVtkm: Encountered unknown DisplaceOperation value!");
    }
}

VISKORES_CONT viskores::cont::DataSet DisplaceFilter::DoExecute(const viskores::cont::DataSet &inputDataset)
{
    auto inputScalar = this->GetFieldFromDataSet(inputDataset);
    auto inputCoords = inputDataset.GetCoordinateSystem(this->GetActiveCoordinateSystemIndex());

    viskores::cont::UnknownArrayHandle outputCoords;

    // TODO: don't hard code dimension...
    this->CastAndCallVecField<3>(inputCoords.GetData(), [&](const auto &coords) {
        using CoordsArrayType = std::decay_t<decltype(coords)>;
        using CoordType = typename CoordsArrayType::ValueType;
        constexpr int N = CoordType::NUM_COMPONENTS;

        this->CastAndCallScalarField(inputScalar.GetData(), [&](const auto &scalars) {
            viskores::cont::ArrayHandle<CoordType> result;
            result.Allocate(coords.GetNumberOfValues());
            viskores::cont::ArrayCopy(coords, result);

            if (m_component == DisplaceComponent::All) {
                for (viskores::IdComponent c = 0; c < N; c++)
                    applyOperation(scalars, result, m_scale, m_operation, c);

                outputCoords = result;
            } else if (m_component == DisplaceComponent::X || m_component == DisplaceComponent::Y ||
                       m_component == DisplaceComponent::Z) {
                viskores::IdComponent c = static_cast<viskores::IdComponent>(m_component);
                if (c >= N)
                    throw viskores::cont::ErrorBadValue(
                        "Error in DisplaceFilter: DisplaceComponent value (" + std::to_string(c) +
                        ") out of bounds for coordinate dimension (" + std::to_string(N) + ")!");
                applyOperation(scalars, result, m_scale, m_operation, c);

                outputCoords = result;

            } else {
                throw viskores::cont::ErrorBadValue(
                    "Error in DisplaceFilter: Encountered unknown DisplaceComponent value!");
            }
        });
    });

    auto outputFieldName = this->GetOutputFieldName();
    if (outputFieldName == "")
        outputFieldName = inputScalar.GetName() + "_displaced";

    return this->CreateResultCoordinateSystem(
        inputDataset, inputDataset.GetCellSet(), inputCoords.GetName(), outputCoords,
        [&](viskores::cont::DataSet &out, const viskores::cont::Field &fieldToPass) {
            out.AddField(viskores::cont::Field(
                fieldToPass.GetName() == this->GetActiveFieldName() ? outputFieldName : fieldToPass.GetName(),
                fieldToPass.GetAssociation(), fieldToPass.GetData()));
        });
}
