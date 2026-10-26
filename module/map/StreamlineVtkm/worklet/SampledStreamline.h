#ifndef VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINE_H
#define VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINE_H

#include <string>
#include <vector>

#include <viskores/Particle.h>
#include <viskores/Types.h>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/DataSet.h>
#include <viskores/cont/Token.h>
#include <viskores/filter/flow/FilterParticleAdvectionSteadyState.h>
#include <viskores/filter/flow/FlowTypes.h>
#include <viskores/filter/flow/worklet/Analysis.h>
#include <viskores/filter/flow/worklet/Field.h>
#include <viskores/filter/flow/worklet/GridEvaluators.h>
#include <viskores/filter/flow/worklet/Termination.h>
#include <viskores/VecVariable.h>

/*
    A variant of viskores::filter::flow::Streamline that samples data fields in its particle
    analysis hooks while accepted integration steps are recorded. Each output partition carries
    values evaluated on the input partition that traced its streamline pieces. Points that could not
    be evaluated there (e.g. at partition boundaries) are NaN.

    Sampled fields are evaluated as viskores::Vec3f internally (missing components are padded with
    zeros, components beyond the third are dropped). Scalar fields are returned as scalars
    (viskores::FloatDefault), all others as viskores::Vec3f.
*/

namespace vistle {

class SampledStreamline;

using SampledFieldType = viskores::worklet::flow::VelocityField<viskores::cont::ArrayHandle<viskores::Vec3f>>;
using SampledEvaluatorType = viskores::worklet::flow::GridEvaluator<SampledFieldType>;

class SampledStreamlineAnalysisExec {
public:
    static constexpr viskores::IdComponent MaxSampledFields = 3;

    using BaseExec = viskores::worklet::flow::StreamlineAnalysisExec<viskores::Particle>;
    using EvaluatorType = viskores::worklet::flow::ExecutionGridEvaluator<SampledFieldType>;
    using ValuesPortal = viskores::cont::ArrayHandle<viskores::Vec3f>::WritePortalType;
    using LengthsPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;
    using ValidityPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;

    VISKORES_CONT SampledStreamlineAnalysisExec(BaseExec base, const std::vector<SampledEvaluatorType> &evaluators,
                                                const std::vector<viskores::cont::ArrayHandle<viskores::Vec3f>> &values,
                                                const viskores::cont::ArrayHandle<viskores::Id> &lengths,
                                                const viskores::cont::ArrayHandle<viskores::Id> &validity,
                                                viskores::IdComponent numFields, viskores::Id maxSteps,
                                                viskores::cont::DeviceAdapterId device, viskores::cont::Token &token);

    VISKORES_EXEC void PreStepAnalyze(viskores::Id index, const viskores::Particle &particle);
    VISKORES_EXEC void Analyze(viskores::Id index, const viskores::Particle &oldParticle,
                               const viskores::Particle &newParticle);

private:
    VISKORES_EXEC void Sample(viskores::Id loc, const viskores::Vec3f &position);

    BaseExec m_base;
    viskores::IdComponent m_numFields;
    viskores::Id m_maxSteps;
    LengthsPortal m_lengths;
    ValidityPortal m_validity;
    // viskores::Vec only supports arithmetic component types, so use plain arrays
    EvaluatorType m_evaluators[MaxSampledFields];
    ValuesPortal m_values[MaxSampledFields];
};

class SampledStreamlineAnalysis: public viskores::worklet::flow::StreamlineAnalysis<viskores::Particle> {
    using Base = viskores::worklet::flow::StreamlineAnalysis<viskores::Particle>;

public:
    using SampledFieldType = vistle::SampledFieldType;
    using EvaluatorType = vistle::SampledEvaluatorType;

    SampledStreamlineAnalysis() = default;
    explicit SampledStreamlineAnalysis(viskores::Id maxSteps);

    // Registers the field with the given name of the dataset for sampling. Returns false (and
    // does nothing) if the field does not exist or has an unsupported association or type.
    bool AddSampledField(const viskores::cont::DataSet &dataset, const std::string &name);

    VISKORES_CONT void UseAsTemplate(const SampledStreamlineAnalysis &other);
    VISKORES_CONT void InitializeAnalysis(const viskores::cont::ArrayHandle<viskores::Particle> &particles);
    VISKORES_CONT SampledStreamlineAnalysisExec PrepareForExecution(viskores::cont::DeviceAdapterId device,
                                                                    viskores::cont::Token &token) const;
    VISKORES_CONT void FinalizeAnalysis(viskores::cont::ArrayHandle<viskores::Particle> &particles);
    static bool MakeDataSet(viskores::cont::DataSet &dataset, const std::vector<SampledStreamlineAnalysis> &results);

private:
    struct Sampler {
        std::string name;
        viskores::IdComponent numComponents = 0;
        EvaluatorType evaluator;
    };

    std::vector<Sampler> m_samplers;
    std::vector<viskores::cont::ArrayHandle<viskores::Vec3f>> m_values;
    viskores::cont::ArrayHandle<viskores::Id> m_streamLengths;
    viskores::cont::ArrayHandle<viskores::Id> m_validity;
    viskores::Id m_maxSteps = 0;
};

} // namespace vistle

namespace viskores {
namespace filter {
namespace flow {

template<>
struct FlowTraits<vistle::SampledStreamline> {
    using ParticleType = viskores::Particle;
    using TerminationType = viskores::worklet::flow::NormalTermination;
    using AnalysisType = vistle::SampledStreamlineAnalysis;
    using ArrayType = viskores::cont::ArrayHandle<viskores::Vec3f>;
    using FieldType = viskores::worklet::flow::VelocityField<ArrayType>;
};

} // namespace flow
} // namespace filter
} // namespace viskores

namespace vistle {

class SampledStreamline: public viskores::filter::flow::FilterParticleAdvectionSteadyState<SampledStreamline> {
public:
    using ParticleType = viskores::filter::flow::FlowTraits<SampledStreamline>::ParticleType;
    using TerminationType = viskores::filter::flow::FlowTraits<SampledStreamline>::TerminationType;
    using AnalysisType = viskores::filter::flow::FlowTraits<SampledStreamline>::AnalysisType;
    using ArrayType = viskores::filter::flow::FlowTraits<SampledStreamline>::ArrayType;
    using FieldType = viskores::filter::flow::FlowTraits<SampledStreamline>::FieldType;

    // Names of up to three point or cell fields to sample during streamline integration.
    void SetSampledFields(const std::vector<std::string> &names) { m_sampledFields = names; }

    // these are called by the filter base class
    FieldType GetField(const viskores::cont::DataSet &data) const;
    TerminationType GetTermination(const viskores::cont::DataSet &data) const;
    AnalysisType GetAnalysis(const viskores::cont::DataSet &data) const;

private:
    std::vector<std::string> m_sampledFields;
};

} // namespace vistle

#endif // VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINE_H
