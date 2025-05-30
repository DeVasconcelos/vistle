#ifndef VISTLE_VTKM_CONVERT_TOPOLOGY_H
#define VISTLE_VTKM_CONVERT_TOPOLOGY_H

#include <vistle/core/object.h>
#include <vistle/core/scalar.h>
#include <vistle/core/triangles.h>
#include <vistle/core/vec.h>

#include <viskores/cont/DataSet.h>

#include "export.h"
#include "convert.h"
#include "module_status.h"

namespace vistle {

Object::ptr vtkmGetTopology(const viskores::cont::DataSet &dataset);
ModuleStatusPtr vtkmSetTopology(viskores::cont::DataSet &vtkmDataset, vistle::Object::const_ptr grid);
} // namespace vistle

#endif
