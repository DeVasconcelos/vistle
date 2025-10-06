#include <viskores/filter/entity_extraction/ExternalFaces.h>

#include "DomainSurfaceVtkm.h"

MODULE_MAIN(DomainSurfaceVtkm)

using namespace vistle;

DomainSurfaceVtkm::DomainSurfaceVtkm(const std::string &name, int moduleID, mpi::communicator comm)
: VtkmModule(name, moduleID, comm, 1, MappedDataHandling::Require)
{}

DomainSurfaceVtkm::~DomainSurfaceVtkm()
{}

std::unique_ptr<viskores::filter::Filter> DomainSurfaceVtkm::setUpFilter() const
{
    return std::make_unique<viskores::filter::entity_extraction::ExternalFaces>();
}
