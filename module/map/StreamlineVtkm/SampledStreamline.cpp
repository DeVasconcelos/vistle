#include <viskores/cont/ErrorFilterExecution.h>
#include <viskores/filter/flow/internal/FilterParticleAdvectionSteadyStateImpl.h>
#include <viskores/filter/flow/worklet/ParticleAdvectionWorklets.h>

#include "SampledStreamline.h"

namespace vistle {

VISKORES_CONT SampledStreamline::FieldType SampledStreamline::GetField(const viskores::cont::DataSet &dataset) const
{
    const auto &fieldNm = this->GetActiveFieldName();
    if (!dataset.HasPointField(fieldNm) && !dataset.HasCellField(fieldNm))
        throw viskores::cont::ErrorFilterExecution("Unsupported field association");
    auto assoc = dataset.GetField(fieldNm).GetAssociation();
    ArrayType arr;
    viskores::cont::ArrayCopyShallowIfPossible(dataset.GetField(fieldNm).GetData(), arr);
    return viskores::worklet::flow::VelocityField<ArrayType>(arr, assoc);
}

VISKORES_CONT SampledStreamline::TerminationType
SampledStreamline::GetTermination(const viskores::cont::DataSet &dataset) const
{
    // dataset not used
    (void)dataset;
    return SampledStreamline::TerminationType(this->NumberOfSteps);
}

VISKORES_CONT SampledStreamline::AnalysisType
SampledStreamline::GetAnalysis(const viskores::cont::DataSet &dataset) const
{
    SampledStreamline::AnalysisType analysis(this->NumberOfSteps);
    for (const auto &name: m_sampledFields)
        if (!analysis.AddSampledField(dataset, name))
            throw viskores::cont::ErrorFilterExecution("Could not sample input field: " + name);

    // by default, sample the active field if no other fields were specified
    if (m_sampledFields.empty() && !analysis.AddSampledField(dataset, this->GetActiveFieldName()))
        throw viskores::cont::ErrorFilterExecution("Could not sample active input field: " +
                                                   this->GetActiveFieldName());
    return analysis;
}

} // namespace vistle

template class viskores::filter::flow::FilterParticleAdvectionSteadyState<vistle::SampledStreamline>;
