#ifndef VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINEANALYSIS_H
#define VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINEANALYSIS_H

#include <string>
#include <vector>

#include <viskores/Particle.h>
#include <viskores/Types.h>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/DataSet.h>
#include <viskores/cont/Token.h>
#include <viskores/filter/flow/worklet/Analysis.h>
#include <viskores/filter/flow/worklet/Field.h>
#include <viskores/filter/flow/worklet/GridEvaluators.h>

namespace vistle {

using SampledFieldType = viskores::worklet::flow::VelocityField<viskores::cont::ArrayHandle<viskores::Vec3f>>;
using SampledEvaluatorType = viskores::worklet::flow::GridEvaluator<SampledFieldType>;

class SampledStreamlineAnalysisExec {
public:
    using BaseExec = viskores::worklet::flow::StreamlineAnalysisExec<viskores::Particle>;
    using EvaluatorType = viskores::worklet::flow::ExecutionGridEvaluator<SampledFieldType>;
    using ValuesPortal = viskores::cont::ArrayHandle<viskores::Vec3f>::WritePortalType;
    struct FieldExec {
        EvaluatorType evaluator;
        ValuesPortal values;
    };
    using FieldsPortal = viskores::cont::ArrayHandle<FieldExec>::ReadPortalType;
    using LengthsPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;
    using ValidityPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;

    VISKORES_CONT SampledStreamlineAnalysisExec(BaseExec base, const viskores::cont::ArrayHandle<FieldExec> &fields,
                                                const viskores::cont::ArrayHandle<viskores::Id> &lengths,
                                                const viskores::cont::ArrayHandle<viskores::Id> &validity,
                                                viskores::Id maxSteps, viskores::cont::DeviceAdapterId device,
                                                viskores::cont::Token &token);

    VISKORES_EXEC void PreStepAnalyze(viskores::Id index, const viskores::Particle &particle);
    VISKORES_EXEC void Analyze(viskores::Id index, const viskores::Particle &oldParticle,
                               const viskores::Particle &newParticle);

private:
    VISKORES_EXEC void Sample(viskores::Id loc, const viskores::Vec3f &position);

    BaseExec m_base;
    viskores::Id m_numFields;
    viskores::Id m_maxSteps;
    LengthsPortal m_lengths;
    ValidityPortal m_validity;
    FieldsPortal m_fields;
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

#endif // VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINEANALYSIS_H
