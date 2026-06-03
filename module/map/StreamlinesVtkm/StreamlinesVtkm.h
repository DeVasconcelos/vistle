#ifndef VISTLE_STREAMLINESVTKM_STREAMLINESVTKM_H
#define VISTLE_STREAMLINESVTKM_STREAMLINESVTKM_H

#include <vistle/vtkm/vtkm_module.h>

class StreamlinesVtkm: public vistle::VtkmModule {
public:
    StreamlinesVtkm(const std::string &name, int moduleID, mpi::communicator comm);
    ~StreamlinesVtkm();

private:
    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;
};

#endif // VISTLE_STREAMLINESVTKM_STREAMLINESVTKM_H