#ifndef VISTLE_STREAMLINEVTKM_SAMPLEDSTREAMLINE_H
#define VISTLE_STREAMLINEVTKM_SAMPLEDSTREAMLINE_H

#include <string>
#include <vector>

#include <viskores/cont/DataSet.h>
#include <viskores/filter/flow/FilterParticleAdvectionSteadyState.h>
#include <viskores/filter/flow/FlowTypes.h>
#include <viskores/filter/flow/worklet/Field.h>
#include <viskores/filter/flow/worklet/Termination.h>

#include "worklet/SampledStreamlineAnalysis.h"

namespace vistle {

class SampledStreamline;

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

/*
    A variant of viskores::filter::flow::Streamline that samples data fields in its particle
    analysis hooks while accepted integration steps are recorded. Each output partition carries
    values evaluated on the input partition that traced its streamline pieces. Points that could not
    be evaluated there (e.g. at partition boundaries) are NaN.

    Sampled fields keep their value type and number of components. Cell fields take the value of the
    cell containing a point, point fields the value of the nearest vertex of that cell (no interpolation).
    Missing values are NaN for floating point fields and 0 for integer fields.
*/
class SampledStreamline: public viskores::filter::flow::FilterParticleAdvectionSteadyState<SampledStreamline> {
public:
    using ParticleType = viskores::filter::flow::FlowTraits<SampledStreamline>::ParticleType;
    using TerminationType = viskores::filter::flow::FlowTraits<SampledStreamline>::TerminationType;
    using AnalysisType = viskores::filter::flow::FlowTraits<SampledStreamline>::AnalysisType;
    using ArrayType = viskores::filter::flow::FlowTraits<SampledStreamline>::ArrayType;
    using FieldType = viskores::filter::flow::FlowTraits<SampledStreamline>::FieldType;

    VISKORES_CONT FieldType GetField(const viskores::cont::DataSet &data) const;

    VISKORES_CONT TerminationType GetTermination(const viskores::cont::DataSet &data) const;

    VISKORES_CONT AnalysisType GetAnalysis(const viskores::cont::DataSet &data) const;

    // Names of point or cell fields to sample during streamline integration.
    void SetSampledFields(const std::vector<std::string> &names) { m_sampledFields = names; }

private:
    std::vector<std::string> m_sampledFields;
};

} // namespace vistle

#endif // VISTLE_STREAMLINEVTKM_SAMPLEDSTREAMLINE_H
