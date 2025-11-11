#ifndef VISTLE_DISPLACEVTKM_FILTER_DISPLACEWORKLET_H
#define VISTLE_DISPLACEVTKM_FILTER_DISPLACEWORKLET_H

#include <viskores/worklet/WorkletMapField.h>

#include <vistle/util/enum.h>

struct BaseDisplaceWorklet {
    viskores::FloatDefault m_scale;

    VISKORES_CONT BaseDisplaceWorklet(): m_scale{1.0f} {};

    VISKORES_CONT BaseDisplaceWorklet(viskores::FloatDefault scale): m_scale{scale} {}
};

struct SetDisplaceWorklet: BaseDisplaceWorklet, viskores::worklet::WorkletMapField {
    VISKORES_CONT SetDisplaceWorklet(): BaseDisplaceWorklet(){};

    VISKORES_CONT SetDisplaceWorklet(viskores::FloatDefault scale): BaseDisplaceWorklet(scale) {}

    using ControlSignature = void(FieldIn scalarField, FieldIn coords, FieldOut result);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const viskores::Vec<T, 3> &scalar, const viskores::Vec<S, 3> &coords,
                                  viskores::Vec<S, 3> &displacedCoord) const
    {
        for (viskores::Id i = 0; i < 3; i++)
            displacedCoord[i] = scalar[i] * this->m_scale;
    }

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const S &coord, S &displacedCoord) const
    {
        displacedCoord = scalar * this->m_scale;
    }
};

struct AddDisplaceWorklet: BaseDisplaceWorklet, viskores::worklet::WorkletMapField {
    VISKORES_CONT AddDisplaceWorklet(): BaseDisplaceWorklet(){};

    VISKORES_CONT AddDisplaceWorklet(viskores::FloatDefault scale): BaseDisplaceWorklet(scale) {}

    using ControlSignature = void(FieldIn scalarField, FieldIn coords, FieldOut result);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const viskores::Vec<T, 3> &scalar, const viskores::Vec<S, 3> &coords,
                                  viskores::Vec<S, 3> &displacedCoord) const
    {
        for (viskores::Id i = 0; i < 3; i++)
            displacedCoord[i] = coords[i] + scalar[i] * this->m_scale;
    }

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const S &coord, S &displacedCoord) const
    {
        displacedCoord = coord + scalar * this->m_scale;
    }
};

struct MultiplyDisplaceWorklet: BaseDisplaceWorklet, viskores::worklet::WorkletMapField {
    VISKORES_CONT MultiplyDisplaceWorklet(): BaseDisplaceWorklet(){};

    VISKORES_CONT MultiplyDisplaceWorklet(viskores::FloatDefault scale): BaseDisplaceWorklet(scale) {}

    using ControlSignature = void(FieldIn scalarField, FieldIn coords, FieldOut result);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const viskores::Vec<T, 3> &scalar, const viskores::Vec<S, 3> &coords,
                                  viskores::Vec<S, 3> &displacedCoord) const
    {
        for (viskores::Id i = 0; i < 3; i++)
            displacedCoord[i] = coords[i] * scalar[i] * this->m_scale;
    }

    template<typename T, typename S>
    VISKORES_EXEC void operator()(const T &scalar, const S &coord, S &displacedCoord) const
    {
        displacedCoord = coord * scalar * this->m_scale;
    }
};

#endif
