#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>

#include <viskores/Math.h>
#include <viskores/TypeList.h>
#include <viskores/cont/Algorithm.h>
#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/ArrayHandleCast.h>
#include <viskores/cont/ArrayHandleCompositeVector.h>
#include <viskores/cont/ArrayHandleConstant.h>
#include <viskores/cont/ArrayHandleExtractComponent.h>
#include <viskores/cont/ErrorFilterExecution.h>
#include <viskores/filter/flow/internal/FilterParticleAdvectionSteadyStateImpl.h>
#include <viskores/filter/flow/worklet/ParticleAdvectionWorklets.h>

#include "SampledStreamline.h"

namespace vistle {

namespace {

// Converts the (flattened) components of an array to a Vec3f array: missing components are padded with zeros,
// components beyond the third are dropped.
bool toVec3f(const viskores::cont::UnknownArrayHandle &in, viskores::cont::ArrayHandle<viskores::Vec3f> &out,
             viskores::IdComponent &numComponents)
{
    numComponents = in.GetNumberOfComponentsFlat();
    const viskores::Id n = in.GetNumberOfValues();

    std::array<viskores::cont::ArrayHandle<viskores::FloatDefault>, 3> components;
    for (auto &c: components)
        viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleConstant<viskores::FloatDefault>(0, n), c);

    bool handled = false;
    auto tryType = [&](auto t) {
        using T = std::decay_t<decltype(t)>;
        if (handled || !in.IsBaseComponentType<T>())
            return;
        handled = true;
        for (viskores::IdComponent c = 0; c < std::min<viskores::IdComponent>(numComponents, 3); ++c) {
            viskores::cont::ArrayCopy(
                viskores::cont::make_ArrayHandleCast<viskores::FloatDefault>(in.ExtractComponent<T>(c)), components[c]);
        }
    };
    viskores::ListForEach(tryType, viskores::TypeListScalarAll{});
    if (!handled)
        return false;

    viskores::cont::ArrayCopy(
        viskores::cont::make_ArrayHandleCompositeVector(components[0], components[1], components[2]), out);
    return true;
}

// unused slots reuse the last valid entry, as the evaluators have no default constructor
SampledStreamlineAnalysisExec::EvaluatorType prepareEvaluator(const std::vector<SampledEvaluatorType> &evaluators,
                                                              std::size_t i, viskores::cont::DeviceAdapterId device,
                                                              viskores::cont::Token &token)
{
    return evaluators[std::min(i, evaluators.size() - 1)].PrepareForExecution(device, token);
}

SampledStreamlineAnalysisExec::ValuesPortal
preparePortal(const std::vector<viskores::cont::ArrayHandle<viskores::Vec3f>> &values, std::size_t i,
              viskores::cont::DeviceAdapterId device, viskores::cont::Token &token)
{
    return values[std::min(i, values.size() - 1)].PrepareForInPlace(device, token);
}

} // namespace

VISKORES_CONT SampledStreamlineAnalysisExec::SampledStreamlineAnalysisExec(
    BaseExec base, const std::vector<SampledEvaluatorType> &evaluators,
    const std::vector<viskores::cont::ArrayHandle<viskores::Vec3f>> &values,
    const viskores::cont::ArrayHandle<viskores::Id> &lengths, const viskores::cont::ArrayHandle<viskores::Id> &validity,
    viskores::IdComponent numFields, viskores::Id maxSteps, viskores::cont::DeviceAdapterId device,
    viskores::cont::Token &token)
: m_base(std::move(base))
, m_numFields(numFields)
, m_maxSteps(maxSteps + 1)
, m_lengths(lengths.PrepareForInPlace(device, token))
, m_validity(validity.PrepareForInPlace(device, token))
, m_evaluators{prepareEvaluator(evaluators, 0, device, token), prepareEvaluator(evaluators, 1, device, token),
               prepareEvaluator(evaluators, 2, device, token)}
, m_values{preparePortal(values, 0, device, token), preparePortal(values, 1, device, token),
           preparePortal(values, 2, device, token)}
{}

SampledStreamlineAnalysis::SampledStreamlineAnalysis(viskores::Id maxSteps): Base(maxSteps), m_maxSteps(maxSteps)
{}

bool SampledStreamlineAnalysis::AddSampledField(const viskores::cont::DataSet &dataset, const std::string &name)
{
    using Association = viskores::cont::Field::Association;

    if (!dataset.HasField(name))
        return false;
    const auto &field = dataset.GetField(name);
    const auto assoc = field.GetAssociation();
    if (assoc != Association::Points && assoc != Association::Cells)
        return false;

    Sampler sampler;
    sampler.name = name;
    viskores::cont::ArrayHandle<viskores::Vec3f> values;
    if (!toVec3f(field.GetData(), values, sampler.numComponents))
        return false;
    sampler.evaluator = EvaluatorType(dataset, SampledFieldType(values, assoc));
    m_samplers.push_back(std::move(sampler));
    return true;
}

VISKORES_CONT void SampledStreamlineAnalysis::UseAsTemplate(const SampledStreamlineAnalysis &other)
{
    Base::UseAsTemplate(other);
    m_samplers = other.m_samplers;
    m_values.clear();
    m_maxSteps = other.m_maxSteps;
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

    m_values.clear();
    for (std::size_t i = 0; i < m_samplers.size(); ++i) {
        viskores::cont::ArrayHandle<viskores::Vec3f> values;
        values.Allocate(numPoints);
        const auto missing = viskores::Vec3f(viskores::Nan<viskores::FloatDefault>());
        viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleConstant<viskores::Vec3f>(missing, numPoints),
                                  values);
        m_values.push_back(values);
    }
}

VISKORES_CONT SampledStreamlineAnalysisExec SampledStreamlineAnalysis::PrepareForExecution(
    viskores::cont::DeviceAdapterId device, viskores::cont::Token &token) const
{
    std::vector<EvaluatorType> evaluators;
    evaluators.reserve(m_samplers.size());
    for (const auto &sampler: m_samplers)
        evaluators.push_back(sampler.evaluator);

    auto base = Base::PrepareForExecution(device, token);
    return SampledStreamlineAnalysisExec(base, evaluators, m_values, m_streamLengths, m_validity,
                                         static_cast<viskores::IdComponent>(m_samplers.size()), m_maxSteps, device,
                                         token);
}

VISKORES_EXEC void SampledStreamlineAnalysisExec::PreStepAnalyze(viskores::Id index, const viskores::Particle &particle)
{
    m_base.PreStepAnalyze(index, particle);
    if (m_lengths.Get(index) == 0) {
        const viskores::Id loc = index * m_maxSteps;
        m_lengths.Set(index, 1);
        m_validity.Set(loc, 1);
        Sample(loc, particle.GetPosition());
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
    Sample(loc, newParticle.GetPosition());
}

VISKORES_EXEC void SampledStreamlineAnalysisExec::Sample(viskores::Id loc, const viskores::Vec3f &position)
{
    for (viskores::Id field = 0; field < m_numFields; ++field) {
        viskores::VecVariable<viskores::Vec3f, 2> value;
        const auto status = m_evaluators[field].Evaluate(position, viskores::FloatDefault(0), value);
        if (status.CheckOk())
            m_values[field].Set(loc, value[0]);
        // Points on or just outside the partition boundary have no local interpolant and stay NaN,
        // so that they can be resolved against the other partitions later.
    }
}

VISKORES_CONT void
SampledStreamlineAnalysis::FinalizeAnalysis(viskores::cont::ArrayHandle<viskores::Particle> &particles)
{
    Base::FinalizeAnalysis(particles);

    for (auto &values: m_values) {
        viskores::cont::ArrayHandle<viskores::Vec3f> compacted;
        viskores::cont::Algorithm::CopyIf(values, m_validity, compacted, Base::IsOne());
        values = compacted;
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
    for (std::size_t s = 0; s < samplers.size(); ++s) {
        viskores::cont::ArrayHandle<viskores::Vec3f> values;
        if (results.size() == 1) {
            values = results.front().m_values[s];
        } else {
            viskores::Id total = 0;
            for (const auto &res: results)
                total += res.m_values[s].GetNumberOfValues();
            values.Allocate(total);
            viskores::Id offset = 0;
            for (const auto &res: results) {
                const auto &part = res.m_values[s];
                viskores::cont::Algorithm::CopySubRange(part, 0, part.GetNumberOfValues(), values, offset);
                offset += part.GetNumberOfValues();
            }
        }

        if (samplers[s].numComponents == 1) {
            viskores::cont::ArrayHandle<viskores::FloatDefault> scalar;
            viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleExtractComponent(values, 0), scalar);
            dataset.AddPointField(samplers[s].name, scalar);
        } else {
            dataset.AddPointField(samplers[s].name, values);
        }
    }
    return true;
}

SampledStreamline::FieldType SampledStreamline::GetField(const viskores::cont::DataSet &dataset) const
{
    const auto &name = this->GetActiveFieldName();
    if (!dataset.HasPointField(name) && !dataset.HasCellField(name))
        throw viskores::cont::ErrorFilterExecution("Unsupported field association");
    auto assoc = dataset.GetField(name).GetAssociation();
    ArrayType array;
    viskores::cont::ArrayCopyShallowIfPossible(dataset.GetField(name).GetData(), array);
    return FieldType(array, assoc);
}

SampledStreamline::TerminationType SampledStreamline::GetTermination(const viskores::cont::DataSet &) const
{
    return TerminationType(this->NumberOfSteps);
}

SampledStreamline::AnalysisType SampledStreamline::GetAnalysis(const viskores::cont::DataSet &dataset) const
{
    if (m_sampledFields.size() > SampledStreamlineAnalysisExec::MaxSampledFields)
        throw viskores::cont::ErrorFilterExecution("At most three input fields can be sampled");

    AnalysisType analysis(this->NumberOfSteps);
    for (const auto &name: m_sampledFields)
        if (!analysis.AddSampledField(dataset, name))
            throw viskores::cont::ErrorFilterExecution("Could not sample input field: " + name);
    // the execution object cannot be default-constructed and needs at least one evaluator
    if (m_sampledFields.empty())
        analysis.AddSampledField(dataset, this->GetActiveFieldName());
    return analysis;
}

} // namespace vistle

template class viskores::filter::flow::FilterParticleAdvectionSteadyState<vistle::SampledStreamline>;
