#include <array>

#include <viskores/TypeList.h>
#include <viskores/VecVariable.h>

#include <viskores/cont/Algorithm.h>
#include <viskores/cont/ArrayCopy.h>
#include <viskores/cont/ArrayHandleCast.h>
#include <viskores/cont/ArrayHandleCompositeVector.h>
#include <viskores/cont/ArrayHandleConstant.h>
#include <viskores/cont/ArrayHandleExtractComponent.h>

#include "SampledStreamlineAnalysis.h"

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

} // namespace

VISKORES_CONT SampledStreamlineAnalysisExec::SampledStreamlineAnalysisExec(
    BaseExec base, const viskores::cont::ArrayHandle<FieldExec> &fields,
    const viskores::cont::ArrayHandle<viskores::Id> &lengths, const viskores::cont::ArrayHandle<viskores::Id> &validity,
    viskores::Id maxSteps, viskores::cont::DeviceAdapterId device, viskores::cont::Token &token)
: m_base(std::move(base))
, m_numFields(fields.GetNumberOfValues())
, m_maxSteps(maxSteps + 1)
, m_lengths(lengths.PrepareForInPlace(device, token))
, m_validity(validity.PrepareForInPlace(device, token))
, m_fields(fields.PrepareForInput(device, token))
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
    std::vector<SampledStreamlineAnalysisExec::FieldExec> fields;
    fields.reserve(m_samplers.size());
    for (std::size_t i = 0; i < m_samplers.size(); ++i) {
        fields.push_back(
            {m_samplers[i].evaluator.PrepareForExecution(device, token), m_values[i].PrepareForInPlace(device, token)});
    }
    auto fieldsHandle = viskores::cont::make_ArrayHandle(fields, viskores::CopyFlag::On);

    auto base = Base::PrepareForExecution(device, token);
    return SampledStreamlineAnalysisExec(base, fieldsHandle, m_streamLengths, m_validity, m_maxSteps, device, token);
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
    for (viskores::Id field = 0; field < m_numFields; field++) {
        const auto sampler = m_fields.Get(field);
        viskores::VecVariable<viskores::Vec3f, 2> value;
        const auto status = sampler.evaluator.Evaluate(position, viskores::FloatDefault(0), value);
        if (status.CheckOk())
            sampler.values.Set(loc, value[0]);
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
    for (std::size_t i = 0; i < samplers.size(); i++) {
        viskores::cont::ArrayHandle<viskores::Vec3f> values;
        if (results.size() == 1) {
            values = results.front().m_values[i];
        } else {
            viskores::Id total = 0;
            for (const auto &res: results)
                total += res.m_values[i].GetNumberOfValues();
            values.Allocate(total);
            viskores::Id offset = 0;
            // concatenate the results from all samplers
            for (const auto &res: results) {
                const auto &part = res.m_values[i];
                viskores::cont::Algorithm::CopySubRange(part, 0, part.GetNumberOfValues(), values, offset);
                offset += part.GetNumberOfValues();
            }
        }

        if (samplers[i].numComponents == 1) {
            viskores::cont::ArrayHandle<viskores::FloatDefault> scalar;
            viskores::cont::ArrayCopy(viskores::cont::make_ArrayHandleExtractComponent(values, 0), scalar);
            dataset.AddPointField(samplers[i].name, scalar);
        } else {
            dataset.AddPointField(samplers[i].name, values);
        }
    }
    return true;
}

} // namespace vistle
