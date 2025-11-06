#include <string>

#include <viskores/cont/ArrayHandleExtractComponent.h>
#include <viskores/cont/ArrayHandleSOA.h>

#include "DisplaceFilter.h"
#include "DisplaceWorklet.h"

VISKORES_CONT DisplaceFilter::DisplaceFilter()
: m_component(DisplaceComponent::X), m_operation(DisplaceOperation::Add), m_scale(1.0f)
{}

VISKORES_CONT viskores::cont::DataSet DisplaceFilter::DoExecute(const viskores::cont::DataSet &inputDataset)
{
    auto inputScalar = this->GetFieldFromDataSet(inputDataset);
    auto inputCoords = inputDataset.GetCoordinateSystem(this->GetActiveCoordinateSystemIndex());

    viskores::cont::UnknownArrayHandle outputCoords;

    this->CastAndCallVecField<3>(inputCoords.GetData(), [&](const auto &coords) {
        using CoordsArrayType = std::decay_t<decltype(coords)>;
        using CoordType = typename CoordsArrayType::ValueType;
        constexpr int N = CoordType::NUM_COMPONENTS;

        this->CastAndCallScalarField(inputScalar.GetData(), [&](const auto &scalars) {
            viskores::cont::ArrayHandle<CoordType> result;

            if (m_component == DisplaceComponent::All) {
                switch (m_operation) {
                case DisplaceOperation::Set:
                    this->Invoke(SetDisplaceWorklet<N>{m_scale}, scalars, coords, result);
                    break;
                case DisplaceOperation::Add:
                    this->Invoke(AddDisplaceWorklet<N>{m_scale}, scalars, coords, result);
                    break;
                case DisplaceOperation::Multiply:
                    this->Invoke(MultiplyDisplaceWorklet<N>{m_scale}, scalars, coords, result);
                    break;
                default:
                    throw viskores::cont::ErrorBadValue(
                        "Error in DisplaceFilter: Encountered unknown DisplaceOperation value!");
                }

                outputCoords = result;
            } else if (m_component == DisplaceComponent::X || m_component == DisplaceComponent::Y ||
                       m_component == DisplaceComponent::Z) {
                viskores::IdComponent c = static_cast<viskores::IdComponent>(m_component);
                if (c >= N)
                    throw viskores::cont::ErrorBadValue(
                        "Error in DisplaceFilter: DisplaceComponent value (" + std::to_string(c) +
                        ") out of bounds for coordinate dimension (" + std::to_string(N) + ")!");

                auto desiredComponent = viskores::cont::ArrayHandleExtractComponent<CoordsArrayType>(coords, c);

                using ComponentType = typename std::decay_t<decltype(desiredComponent)>::ValueType;
                viskores::cont::ArrayHandle<ComponentType> result;

                switch (m_operation) {
                case DisplaceOperation::Set:
                    this->Invoke(SetDisplaceWorklet<1>{m_scale}, scalars, desiredComponent, result);
                    break;
                case DisplaceOperation::Add:
                    this->Invoke(AddDisplaceWorklet<1>{m_scale}, scalars, desiredComponent, result);
                    break;
                case DisplaceOperation::Multiply:
                    this->Invoke(MultiplyDisplaceWorklet<1>{m_scale}, scalars, desiredComponent, result);
                    break;
                default:
                    throw viskores::cont::ErrorBadValue(
                        "Error in DisplaceFilter: Encountered unknown DisplaceOperation value!");
                }

                viskores::cont::ArrayHandle<CoordType> resultCoords;
                resultCoords.Allocate(coords.GetNumberOfValues());
                auto resultCoordsPortal = resultCoords.WritePortal();
                auto coordsPortal = coords.ReadPortal();
                auto resultPortal = result.ReadPortal();

                for (auto i = 0; i < coords.GetNumberOfValues(); i++) {
                    if (m_component == DisplaceComponent::X) {
                        resultCoordsPortal.Set(i, viskores::Vec<ComponentType, N>{resultPortal.Get(i),
                                                                                  coordsPortal.Get(i)[1],
                                                                                  coordsPortal.Get(i)[2]});
                    } else if (m_component == DisplaceComponent::Y) {
                        resultCoordsPortal.Set(i, viskores::Vec<ComponentType, N>{coordsPortal.Get(i)[0],
                                                                                  resultPortal.Get(i),
                                                                                  coordsPortal.Get(i)[2]});
                    } else if (m_component == DisplaceComponent::Z) {
                        resultCoordsPortal.Set(i, viskores::Vec<ComponentType, N>{coordsPortal.Get(i)[0],
                                                                                  coordsPortal.Get(i)[1],
                                                                                  resultPortal.Get(i)});
                    }
                }

                outputCoords = resultCoords;

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
