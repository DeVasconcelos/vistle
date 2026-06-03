#include <viskores/filter/flow/Streamline.h>

#include "StreamlinesVtkm.h"

MODULE_MAIN(StreamlinesVtkm)

using namespace vistle;

StreamlinesVtkm::StreamlinesVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 3, MappedDataHandling::Require)
{}

StreamlinesVtkm::~StreamlinesVtkm()
{}

std::unique_ptr<viskores::filter::Filter> StreamlinesVtkm::setUpFilter() const
{
    auto filter = std::make_unique<viskores::filter::flow::Streamline>();

    viskores::cont::ArrayHandle<viskores::Particle> seedArray;
    seedArray.Allocate(2);
    seedArray.WritePortal().Set(0, viskores::Particle({0, 0, 0}, 0));
    seedArray.WritePortal().Set(1, viskores::Particle({1, 1, 1}, 1));

    filter->SetStepSize(0.1f);
    filter->SetNumberOfSteps(100);
    filter->SetSeeds(seedArray);

    return filter;
}
