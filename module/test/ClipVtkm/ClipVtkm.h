#ifndef VISTLE_VTKM_CLIP_H
#define VISTLE_VTKM_CLIP_H

#include <vistle/vtkm/ImplicitFunctionController.h>
#include <vistle/vtkm/VtkmModule.h>

class ClipVtkm: public VtkmModule {
public:
    ClipVtkm(const std::string &name, int moduleID, mpi::communicator comm);
    ~ClipVtkm();

private:
    vistle::IntParameter *m_flip = nullptr;

    void runFilter(vtkm::cont::DataSet &input, vtkm::cont::DataSet &output) const override;
    bool changeParameter(const vistle::Parameter *param) override;

    ImplicitFunctionController m_implFuncControl;
};

#endif // VISTLE_VTKM_ISOSURFACE_H
