set(HEADERS worklet/GenerateSeeds.h worklet/SampledStreamlineAnalysis.h SampledStreamline.h StreamlineVtkm.h)
set(SOURCES worklet/GenerateSeeds.cpp worklet/SampledStreamlineAnalysis.cpp SampledStreamline.cpp StreamlineVtkm.cpp)

add_vtkm_module(StreamlineVtkm "create streamlines for a vector field" ${HEADERS} DEVICE_SOURCES ${SOURCES})
