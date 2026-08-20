#include "GenerateSeeds.h"

void GenerateSeedsOnPlaneWorklet::operator()(viskores::Particle &seed, const viskores::Id &index) const
{
    viskores::Id i = index / m_n1;
    viskores::Id j = index % m_n1;

    viskores::FloatDefault s0 = viskores::FloatDefault(1) / (m_n0 - 1);
    viskores::FloatDefault s1 = viskores::FloatDefault(1) / (m_n1 - 1);
    auto point = m_startpoint1 + m_parallelSpan * s0 * i + m_orthogonalSpan * s1 * j;
    seed = viskores::Particle({point[0], point[1], point[2]}, index);
}

void GenerateSeedsOnLineWorklet::operator()(viskores::Particle &seed, const viskores::Id &index) const
{
    seed = viskores::Particle({m_startpoint1[0] + index * m_delta[0], m_startpoint1[1] + index * m_delta[1],
                               m_startpoint1[2] + index * m_delta[2]},
                              index);
}