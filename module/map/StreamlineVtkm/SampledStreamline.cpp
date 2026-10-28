#include <viskores/cont/ErrorFilterExecution.h>
#include <viskores/filter/flow/internal/FilterParticleAdvectionSteadyStateImpl.h>
#include <viskores/filter/flow/worklet/ParticleAdvectionWorklets.h>

#include "SampledStreamline.h"

namespace vistle {

// same as streamline
SampledStreamline::FieldType SampledStreamline::GetField(const viskores::cont::DataSet &dataset) const
{
    const auto &fieldNm = this->GetActiveFieldName();
    if (!dataset.HasPointField(fieldNm) && !dataset.HasCellField(fieldNm))
        throw viskores::cont::ErrorFilterExecution("Unsupported field association");
    auto assoc = dataset.GetField(fieldNm).GetAssociation();
    ArrayType arr;
    viskores::cont::ArrayCopyShallowIfPossible(dataset.GetField(fieldNm).GetData(), arr);
    return viskores::worklet::flow::VelocityField<ArrayType>(arr, assoc);
}

// also the same as streamline
SampledStreamline::TerminationType SampledStreamline::GetTermination(const viskores::cont::DataSet &dataset) const
{
    // dataset not used
    (void)dataset;
    return SampledStreamline::TerminationType(this->NumberOfSteps);
}

// different from streamline: adds the sampled fields to the analysis object
SampledStreamline::AnalysisType SampledStreamline::GetAnalysis(const viskores::cont::DataSet &dataset) const
{
    if (m_sampledFields.size() > SampledStreamlineAnalysisExec::MaxSampledFields)
        throw viskores::cont::ErrorFilterExecution("At most three input fields can be sampled");

    SampledStreamline::AnalysisType analysis(this->NumberOfSteps);
    for (const auto &name: m_sampledFields)
        if (!analysis.AddSampledField(dataset, name))
            throw viskores::cont::ErrorFilterExecution("Could not sample input field: " + name);
    // by default, sample the active field if no other fields were specified
    if (m_sampledFields.empty())
        analysis.AddSampledField(dataset, this->GetActiveFieldName());
    return analysis;
}

} // namespace vistle

template class viskores::filter::flow::FilterParticleAdvectionSteadyState<vistle::SampledStreamline>;
