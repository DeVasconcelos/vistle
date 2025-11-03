#ifndef VISTLE_DISPLACEVTKM_FILTER_DISPLACEWORKLET_H
#define VISTLE_DISPLACEVTKM_FILTER_DISPLACEWORKLET_H

#include <viskores/worklet/WorkletMapField.h>

#include <vistle/util/enum.h>

template<viskores::IdComponent N>
struct BaseDisplaceWorklet {
    viskores::FloatDefault m_scale;

    VISKORES_CONT BaseDisplaceWorklet(): m_scale{1.0f} {};

    VISKORES_CONT BaseDisplaceWorklet(viskores::FloatDefault scale): m_scale{scale} {}
};

template<viskores::IdComponent N>
struct SetDisplaceWorklet: BaseDisplaceWorklet<N>, viskores::worklet::WorkletMapField {
    VISKORES_CONT SetDisplaceWorklet(): BaseDisplaceWorklet<N>(){};

    VISKORES_CONT SetDisplaceWorklet(viskores::FloatDefault scale): BaseDisplaceWorklet<N>(scale) {}

    using ControlSignature = void(FieldIn scalarField, FieldIn coords, FieldOut result);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const viskores::Vec<S, N> &coord,
                                  viskores::Vec<S, N> &displacedCoord) const
    {
        for (auto c = 0; c < N; c++)
            displacedCoord[c] = scalar * this->m_scale;
    }

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const S &coord, S &displacedCoord) const
    {
        displacedCoord = scalar * this->m_scale;
    }
};

template<viskores::IdComponent N>
struct AddDisplaceWorklet: BaseDisplaceWorklet<N>, viskores::worklet::WorkletMapField {
    VISKORES_CONT AddDisplaceWorklet(): BaseDisplaceWorklet<N>(){};

    VISKORES_CONT AddDisplaceWorklet(viskores::FloatDefault scale): BaseDisplaceWorklet<N>(scale) {}

    using ControlSignature = void(FieldIn scalarField, FieldIn coords, FieldOut result);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const viskores::Vec<S, N> &coord,
                                  viskores::Vec<S, N> &displacedCoord) const
    {
        for (auto c = 0; c < N; c++)
            displacedCoord[c] = coord[c] + scalar * this->m_scale;
    }

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const S &coord, S &displacedCoord) const
    {
        displacedCoord = coord + scalar * this->m_scale;
    }
};

template<viskores::IdComponent N>
struct MultiplyDisplaceWorklet: BaseDisplaceWorklet<N>, viskores::worklet::WorkletMapField {
    VISKORES_CONT MultiplyDisplaceWorklet(): BaseDisplaceWorklet<N>(){};

    VISKORES_CONT MultiplyDisplaceWorklet(viskores::FloatDefault scale): BaseDisplaceWorklet<N>(scale) {}

    using ControlSignature = void(FieldIn scalarField, FieldIn coords, FieldOut result);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const viskores::Vec<S, N> &coord,
                                  viskores::Vec<S, N> &displacedCoord) const
    {
        for (auto c = 0; c < N; c++)
            displacedCoord[c] = coord[c] * scalar * this->m_scale;
    }

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const S &coord, S &displacedCoord) const
    {
        displacedCoord = coord * scalar * this->m_scale;
    }
};

#endif
