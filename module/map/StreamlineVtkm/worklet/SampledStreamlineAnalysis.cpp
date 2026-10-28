#include <limits>
#include <type_traits>

#include <viskores/CellClassification.h>
#include <viskores/TypeList.h>
#include <viskores/VecVariable.h>
#include <viskores/exec/ParametricCoordinates.h>

#include <viskores/cont/Algorithm.h>
#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/ArrayHandleConstant.h>
#include <viskores/cont/Invoker.h>
#include <viskores/exec/CellInterpolate.h>
#include <viskores/worklet/WorkletMapField.h>

#include "SampledStreamlineAnalysis.h"

namespace vistle {
namespace {

// Finds the vertex of the cell containing a sampling position that is nearest to it in parametric space
struct NearestVertexWorklet: viskores::worklet::WorkletMapField {
    using ControlSignature = void(FieldIn cellId, FieldIn parametric, ExecObject helper, FieldOut pointId);
    using ExecutionSignature = void(_1, _2, _3, _4);

    VISKORES_EXEC void operator()(viskores::Id cellId, const viskores::Vec3f &parametric,
                                  const viskores::exec::CellInterpolationHelper &helper, viskores::Id &pointId) const
    {
        pointId = -1;
        if (cellId < 0)
            return;

        viskores::UInt8 shape;
        viskores::IdComponent numVerts;
        viskores::VecVariable<viskores::Id, 8> pointIds;
        helper.GetCellInfo(cellId, shape, numVerts, pointIds);

        viskores::FloatDefault best = std::numeric_limits<viskores::FloatDefault>::max();
        for (viskores::IdComponent v = 0; v < numVerts; ++v) {
            viskores::Vec3f pcoords;
            if (viskores::exec::ParametricCoordinatesPoint(numVerts, v, viskores::CellShapeTagGeneric(shape),
                                                           pcoords) != viskores::ErrorCode::Success)
                continue;
            const auto dist = viskores::MagnitudeSquared(pcoords - parametric);
            if (dist < best) {
                best = dist;
                pointId = pointIds[v];
            }
        }
    }
};

struct InterpolatePointFieldWorklet: viskores::worklet::WorkletMapField {
    using ControlSignature = void(FieldIn cellId, FieldIn parametric, ExecObject helper, WholeArrayIn source,
                                  FieldOut value);
    using ExecutionSignature = void(_1, _2, _3, _4, _5);

    template<typename SourcePortal, typename T>
    VISKORES_EXEC void operator()(viskores::Id cellId, const viskores::Vec3f &parametric,
                                  const viskores::exec::CellInterpolationHelper &helper, const SourcePortal &source,
                                  T &value) const
    {
        value = viskores::Nan<T>();
        if (cellId < 0)
            return;

        viskores::UInt8 shape;
        viskores::IdComponent numVerts;
        viskores::VecVariable<viskores::Id, 8> pointIds;
        helper.GetCellInfo(cellId, shape, numVerts, pointIds);
        if (viskores::exec::CellInterpolate(pointIds, source, parametric, viskores::CellShapeTagGeneric(shape),
                                            value) != viskores::ErrorCode::Success)
            value = viskores::Nan<T>();
    }
};

// Copies the entries of source selected by index to target; the value for an index of -1 is NaN (or 0 for
// integer types)
struct GatherWorklet: viskores::worklet::WorkletMapField {
    using ControlSignature = void(FieldIn index, WholeArrayIn source, FieldInOut target);
    using ExecutionSignature = void(_1, _2, _3);

    template<typename SourcePortal, typename T>
    VISKORES_EXEC void operator()(viskores::Id index, const SourcePortal &source, T &target) const
    {
        if (index >= 0) {
            target = source.Get(index);
        } else if constexpr (std::is_floating_point<T>::value) {
            target = viskores::Nan<T>();
        } else {
            target = T(0);
        }
    }
};

// Gathers the values of data selected by index (-1: missing) into a new array of the same type as data
viskores::cont::UnknownArrayHandle gather(const viskores::cont::UnknownArrayHandle &data,
                                          const viskores::cont::ArrayHandle<viskores::Id> &indices)
{
    auto result = data.NewInstanceBasic();
    result.Allocate(indices.GetNumberOfValues());

    bool handled = false;
    auto tryType = [&](auto t) {
        using T = std::decay_t<decltype(t)>;
        if (handled || !data.IsBaseComponentType<T>())
            return;
        handled = true;
        viskores::cont::Invoker invoke;
        for (viskores::IdComponent c = 0; c < data.GetNumberOfComponentsFlat(); ++c) {
            // extracting without copying: the stride array aliases the buffer of result, so it is updated in place
            auto target = result.ExtractComponent<T>(c, viskores::CopyFlag::Off);
            invoke(GatherWorklet{}, indices, data.ExtractComponent<T>(c), target);
        }
    };
    viskores::ListForEach(tryType, viskores::TypeListScalarAll{});
    return result;
}

viskores::cont::UnknownArrayHandle interpolatePointField(const viskores::cont::UnknownArrayHandle &data,
                                                         const viskores::cont::ArrayHandle<viskores::Id> &cellIds,
                                                         const viskores::cont::ArrayHandle<viskores::Vec3f> &parametric,
                                                         const viskores::cont::ArrayHandle<viskores::Id> &validity,
                                                         const viskores::cont::CellInterpolationHelper &helper)
{
    auto result = data.NewInstanceBasic();
    result.Allocate(cellIds.GetNumberOfValues());
    auto compacted = data.NewInstanceBasic();
    compacted.Allocate(validity.GetNumberOfValues());

    bool handled = false;
    bool interpolated = false;
    auto tryType = [&](auto t) {
        using T = std::decay_t<decltype(t)>;
        if (handled || !data.IsBaseComponentType<T>())
            return;
        handled = true;
        if constexpr (std::is_floating_point<T>::value) {
            interpolated = true;
            viskores::cont::Invoker invoke;
            for (viskores::IdComponent c = 0; c < data.GetNumberOfComponentsFlat(); ++c) {
                auto target = result.ExtractComponent<T>(c, viskores::CopyFlag::Off);
                invoke(InterpolatePointFieldWorklet{}, cellIds, parametric, helper, data.ExtractComponent<T>(c),
                       target);
                auto inputComponent = result.ExtractComponent<T>(c);
                viskores::cont::ArrayHandle<T> compactedComponent;
                viskores::cont::Algorithm::CopyIf(inputComponent, validity, compactedComponent);
                auto outputComponent = compacted.ExtractComponent<T>(c, viskores::CopyFlag::Off);
                viskores::cont::ArrayCopy(compactedComponent, outputComponent);
            }
        }
    };
    viskores::ListForEach(tryType, viskores::TypeListScalarAll{});
    if (!handled || !interpolated)
        return {};

    return compacted;
}

} // namespace

VISKORES_CONT SampledStreamlineAnalysisExec::SampledStreamlineAnalysisExec(
    BaseExec base, const viskores::cont::CellLocatorGeneral &locator,
    const viskores::cont::ArrayHandle<viskores::UInt8> &ghostCells,
    const viskores::cont::ArrayHandle<viskores::Id> &cellIds,
    const viskores::cont::ArrayHandle<viskores::Vec3f> &parametric,
    const viskores::cont::ArrayHandle<viskores::Id> &lengths, const viskores::cont::ArrayHandle<viskores::Id> &validity,
    viskores::Id maxSteps, bool sample, viskores::cont::DeviceAdapterId device, viskores::cont::Token &token)
: m_base(std::move(base))
, m_sample(sample)
, m_haveGhostCells(ghostCells.GetNumberOfValues() > 0)
, m_maxSteps(maxSteps + 1)
, m_locator(locator.PrepareForExecution(device, token))
, m_ghostCells(ghostCells.PrepareForInput(device, token))
, m_cellIds(cellIds.PrepareForInPlace(device, token))
, m_parametric(parametric.PrepareForInPlace(device, token))
, m_lengths(lengths.PrepareForInPlace(device, token))
, m_validity(validity.PrepareForInPlace(device, token))
{}

SampledStreamlineAnalysis::SampledStreamlineAnalysis(viskores::Id maxSteps): Base(maxSteps), m_maxSteps(maxSteps)
{}

void SampledStreamlineAnalysis::SetUpLocator(const viskores::cont::DataSet &dataset)
{
    m_locator.SetCoordinates(dataset.GetCoordinateSystem());
    m_locator.SetCellSet(dataset.GetCellSet());
    m_locator.Update();
    m_interpolationHelper = viskores::cont::CellInterpolationHelper(dataset.GetCellSet());
    m_ghostCells = viskores::cont::ArrayHandle<viskores::UInt8>();
    if (dataset.HasGhostCellField())
        viskores::cont::ArrayCopyShallowIfPossible(dataset.GetGhostCellField().GetData(), m_ghostCells);
    m_haveLocator = true;
}

bool SampledStreamlineAnalysis::AddSampledField(const viskores::cont::DataSet &dataset, const std::string &name)
{
    using Association = viskores::cont::Field::Association;

    if (!dataset.HasField(name))
        return false;
    const auto &field = dataset.GetField(name);
    const auto assoc = field.GetAssociation();
    if (assoc != Association::Points && assoc != Association::Cells)
        return false;

    bool supported = false;
    viskores::ListForEach(
        [&](auto t) { supported = supported || field.GetData().IsBaseComponentType<std::decay_t<decltype(t)>>(); },
        viskores::TypeListScalarAll{});
    if (!supported)
        return false;

    if (!m_haveLocator)
        SetUpLocator(dataset);

    Sampler sampler;
    sampler.name = name;
    sampler.association = assoc;
    sampler.data = field.GetData();
    m_samplers.push_back(std::move(sampler));
    return true;
}

VISKORES_CONT void SampledStreamlineAnalysis::UseAsTemplate(const SampledStreamlineAnalysis &other)
{
    Base::UseAsTemplate(other);
    m_samplers = other.m_samplers;
    m_values.clear();
    m_maxSteps = other.m_maxSteps;
    m_haveLocator = other.m_haveLocator;
    m_locator = other.m_locator;
    m_interpolationHelper = other.m_interpolationHelper;
    m_ghostCells = other.m_ghostCells;
}

VISKORES_CONT void
SampledStreamlineAnalysis::InitializeAnalysis(const viskores::cont::ArrayHandle<viskores::Particle> &particles)
{
    Base::InitializeAnalysis(particles);

    const viskores::Id numParticles = particles.GetNumberOfValues();
    const viskores::Id numPoints = numParticles * (m_maxSteps + 1);
    m_streamLengths.Allocate(numParticles);
    viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleConstant<viskores::Id>(0, numParticles), m_streamLengths);
    m_validity.Allocate(numPoints);
    viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleConstant<viskores::Id>(0, numPoints), m_validity);
    m_cellIds.Allocate(numPoints);
    viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleConstant<viskores::Id>(-1, numPoints), m_cellIds);
    m_parametric.Allocate(numPoints);
    viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleConstant<viskores::Vec3f>(viskores::Vec3f(0), numPoints),
                              m_parametric);
    m_values.clear();
}

VISKORES_CONT SampledStreamlineAnalysisExec SampledStreamlineAnalysis::PrepareForExecution(
    viskores::cont::DeviceAdapterId device, viskores::cont::Token &token) const
{
    auto base = Base::PrepareForExecution(device, token);
    return SampledStreamlineAnalysisExec(base, m_locator, m_ghostCells, m_cellIds, m_parametric, m_streamLengths,
                                         m_validity, m_maxSteps, m_haveLocator && !m_samplers.empty(), device, token);
}

VISKORES_EXEC void SampledStreamlineAnalysisExec::PreStepAnalyze(viskores::Id index, const viskores::Particle &particle)
{
    m_base.PreStepAnalyze(index, particle);
    if (m_lengths.Get(index) == 0) {
        const viskores::Id loc = index * m_maxSteps;
        m_lengths.Set(index, 1);
        m_validity.Set(loc, 1);
        Locate(loc, particle.GetPosition());
    }
}

VISKORES_EXEC void SampledStreamlineAnalysisExec::Analyze(viskores::Id index, const viskores::Particle &oldParticle,
                                                          const viskores::Particle &newParticle)
{
    m_base.Analyze(index, oldParticle, newParticle);
    const viskores::Id length = m_lengths.Get(index);
    const viskores::Id loc = index * m_maxSteps + length;
    m_lengths.Set(index, length + 1);
    m_validity.Set(loc, 1);
    Locate(loc, newParticle.GetPosition());
}

VISKORES_EXEC void SampledStreamlineAnalysisExec::Locate(viskores::Id loc, const viskores::Vec3f &position)
{
    if (!m_sample)
        return;
    viskores::Id cellId = -1;
    viskores::Vec3f parametric;
    m_locator.FindCell(position, cellId, parametric);
    if (cellId != -1 && m_haveGhostCells && m_ghostCells.Get(cellId) == viskores::CellClassification::Ghost)
        cellId = -1;
    m_cellIds.Set(loc, cellId);
    m_parametric.Set(loc, parametric);
}

VISKORES_CONT void
SampledStreamlineAnalysis::FinalizeAnalysis(viskores::cont::ArrayHandle<viskores::Particle> &particles)
{
    Base::FinalizeAnalysis(particles);

    viskores::cont::ArrayHandle<viskores::Id> pointIds;
    bool haveNearestVertices = false;
    viskores::cont::ArrayHandle<viskores::Id> cellIds;
    viskores::cont::Algorithm::CopyIf(m_cellIds, m_validity, cellIds, Base::IsOne());

    m_values.clear();
    for (const auto &sampler: m_samplers) {
        if (sampler.association == viskores::cont::Field::Association::Cells) {
            m_values.push_back(gather(sampler.data, cellIds));
            continue;
        }
        if (sampler.data.IsBaseComponentType<viskores::Float32>() ||
            sampler.data.IsBaseComponentType<viskores::Float64>()) {
            m_values.push_back(
                interpolatePointField(sampler.data, m_cellIds, m_parametric, m_validity, m_interpolationHelper));
            continue;
        }
        if (!haveNearestVertices) {
            viskores::cont::Invoker invoke;
            viskores::cont::ArrayHandle<viskores::Id> allPointIds;
            invoke(NearestVertexWorklet{}, m_cellIds, m_parametric, m_interpolationHelper, allPointIds);
            viskores::cont::Algorithm::CopyIf(allPointIds, m_validity, pointIds, Base::IsOne());
            haveNearestVertices = true;
        }
        m_values.push_back(gather(sampler.data, pointIds));
    }
}

bool SampledStreamlineAnalysis::MakeDataSet(viskores::cont::DataSet &dataset,
                                            const std::vector<SampledStreamlineAnalysis> &results)
{
    if (results.empty())
        return false;

    // the base class concatenates the points of all results in order, so we do the same with the values
    std::vector<Base> bases(results.begin(), results.end());
    if (!Base::MakeDataSet(dataset, bases))
        return false;

    const auto &samplers = results.front().m_samplers;
    for (std::size_t i = 0; i < samplers.size(); i++) {
        viskores::cont::UnknownArrayHandle values = results.front().m_values[i];
        if (results.size() > 1) {
            viskores::Id total = 0;
            for (const auto &res: results)
                total += res.m_values[i].GetNumberOfValues();
            values = results.front().m_values[i].NewInstanceBasic();
            values.Allocate(total);

            bool handled = false;
            auto tryType = [&](auto t) {
                using T = std::decay_t<decltype(t)>;
                if (handled || !values.IsBaseComponentType<T>())
                    return;
                handled = true;
                for (viskores::IdComponent c = 0; c < values.GetNumberOfComponentsFlat(); ++c) {
                    auto out = values.ExtractComponent<T>(c, viskores::CopyFlag::Off).WritePortal();
                    viskores::Id offset = 0;
                    for (const auto &res: results) {
                        auto in = res.m_values[i].ExtractComponent<T>(c).ReadPortal();
                        for (viskores::Id k = 0; k < in.GetNumberOfValues(); ++k)
                            out.Set(offset + k, in.Get(k));
                        offset += in.GetNumberOfValues();
                    }
                }
            };
            viskores::ListForEach(tryType, viskores::TypeListScalarAll{});
            if (!handled)
                return false;
        }

        dataset.AddField(viskores::cont::Field(samplers[i].name, viskores::cont::Field::Association::Points, values));
    }
    return true;
}

} // namespace vistle
