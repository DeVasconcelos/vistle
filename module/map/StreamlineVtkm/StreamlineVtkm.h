#ifndef VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
#define VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H

#include <vistle/vtkm/vtkm_module.h>

class StreamlineVtkm: public vistle::VtkmModule {
public:
    StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm);

private:
    vistle::IntParameter *m_numberOfSteps;
    vistle::FloatParameter *m_stepSize;

    vistle::VectorParameter *m_startPoint1;
    vistle::VectorParameter *m_startPoint2;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const override;
};

#endif // VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
