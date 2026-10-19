set(HEADERS worklet/GenerateSeeds.h worklet/SampledStreamline.h StreamlineVtkm.h)
set(SOURCES worklet/GenerateSeeds.cpp worklet/SampledStreamline.cpp StreamlineVtkm.cpp)

add_vtkm_module(StreamlineVtkm "create streamlines for a vector field" ${HEADERS} DEVICE_SOURCES ${SOURCES})
