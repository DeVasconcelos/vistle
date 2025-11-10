set(HEADERS DisplaceWorklet.h ../MapDrape/MapDrape.h)
set(SOURCES ../MapDrape/MapDrape.cpp)

add_vtkm_module(DisplaceVtkm "move vertices according to mapped data of first input" ${HEADERS} ${SOURCES})
target_compile_definitions(DisplaceVtkm PRIVATE DISPLACEVTKM)

if(VISTLE_USE_CUDA AND VISTLE_USE_SEPARABLE_COMPILATION)
    viskores_add_target_information(vistle_viskores MODIFY_CUDA_FLAGS DEVICE_SOURCES ${SOURCES})
endif()
