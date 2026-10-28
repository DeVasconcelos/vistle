#ifndef VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINEANALYSIS_H
#define VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINEANALYSIS_H

#include <string>
#include <vector>

#include <viskores/Particle.h>
#include <viskores/Types.h>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/CellLocatorGeneral.h>
#include <viskores/cont/DataSet.h>
#include <viskores/cont/Token.h>
#include <viskores/cont/UnknownArrayHandle.h>
#include <viskores/filter/flow/worklet/Analysis.h>
#include <viskores/filter/flow/worklet/CellInterpolationHelper.h>

namespace vistle {

// Records where each accepted step lies in the mesh (cell id and parametric coordinates). The field values are
// gathered afterwards in their original type, so that no conversion to float is necessary.
class SampledStreamlineAnalysisExec {
public:
    using BaseExec = viskores::worklet::flow::StreamlineAnalysisExec<viskores::Particle>;
    using LocatorExec = viskores::cont::CellLocatorGeneral::ExecObjType;
    using CellIdPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;
    using ParametricPortal = viskores::cont::ArrayHandle<viskores::Vec3f>::WritePortalType;
    using LengthsPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;
    using ValidityPortal = viskores::cont::ArrayHandle<viskores::Id>::WritePortalType;
    using GhostPortal = viskores::cont::ArrayHandle<viskores::UInt8>::ReadPortalType;

    VISKORES_CONT SampledStreamlineAnalysisExec(BaseExec base, const viskores::cont::CellLocatorGeneral &locator,
                                                const viskores::cont::ArrayHandle<viskores::UInt8> &ghostCells,
                                                const viskores::cont::ArrayHandle<viskores::Id> &cellIds,
                                                const viskores::cont::ArrayHandle<viskores::Vec3f> &parametric,
                                                const viskores::cont::ArrayHandle<viskores::Id> &lengths,
                                                const viskores::cont::ArrayHandle<viskores::Id> &validity,
                                                viskores::Id maxSteps, bool sample,
                                                viskores::cont::DeviceAdapterId device, viskores::cont::Token &token);

    VISKORES_EXEC void PreStepAnalyze(viskores::Id index, const viskores::Particle &particle);
    VISKORES_EXEC void Analyze(viskores::Id index, const viskores::Particle &oldParticle,
                               const viskores::Particle &newParticle);

private:
    VISKORES_EXEC void Locate(viskores::Id loc, const viskores::Vec3f &position);

    BaseExec m_base;
    bool m_sample;
    bool m_haveGhostCells;
    viskores::Id m_maxSteps;
    LocatorExec m_locator;
    GhostPortal m_ghostCells;
    CellIdPortal m_cellIds;
    ParametricPortal m_parametric;
    LengthsPortal m_lengths;
    ValidityPortal m_validity;
};

class SampledStreamlineAnalysis: public viskores::worklet::flow::StreamlineAnalysis<viskores::Particle> {
    using Base = viskores::worklet::flow::StreamlineAnalysis<viskores::Particle>;

public:
    SampledStreamlineAnalysis() = default;
    explicit SampledStreamlineAnalysis(viskores::Id maxSteps);

    // Registers the field with the given name of the dataset for sampling. Returns false (and
    // does nothing) if the field does not exist or has an unsupported association or type.
    // Values are sampled in the type of the field: cell fields take the value of the containing cell, point
    // fields the value of the vertex of the containing cell that is nearest to the sampling position.
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
        viskores::cont::Field::Association association = viskores::cont::Field::Association::Points;
        viskores::cont::UnknownArrayHandle data;
    };

    void SetUpLocator(const viskores::cont::DataSet &dataset);

    std::vector<Sampler> m_samplers;
    std::vector<viskores::cont::UnknownArrayHandle> m_values;
    bool m_haveLocator = false;
    viskores::cont::CellLocatorGeneral m_locator;
    viskores::cont::CellInterpolationHelper m_interpolationHelper;
    viskores::cont::ArrayHandle<viskores::UInt8> m_ghostCells;
    viskores::cont::ArrayHandle<viskores::Id> m_cellIds;
    viskores::cont::ArrayHandle<viskores::Vec3f> m_parametric;
    viskores::cont::ArrayHandle<viskores::Id> m_streamLengths;
    viskores::cont::ArrayHandle<viskores::Id> m_validity;
    viskores::Id m_maxSteps = 0;
};

} // namespace vistle

#endif // VISTLE_STREAMLINEVTKM_WORKLET_SAMPLEDSTREAMLINEANALYSIS_H
