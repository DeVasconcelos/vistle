#ifndef COORDS_H
#define COORDS_H

#include "scalar.h"
#include "shm.h"
#include "object.h"
#include "vec.h"
#include "normals.h"
#include "export.h"

namespace vistle {

class V_COREEXPORT Coords: public Vec<Scalar,3> {
   V_OBJECT(Coords);

 public:
   typedef Vec<Scalar,3> Base;

   Coords(const Index numVertices,
             const Meta &meta=Meta());

   Index getNumCoords() const;
   Index getNumVertices() const;
   Normals::const_ptr normals() const;
   void setNormals(Normals::const_ptr normals);

   V_DATA_BEGIN(Coords);
      shm_obj_ref<Normals> normals;

      Data(const Index numVertices = 0,
            Type id = UNKNOWN, const std::string & name = "",
            const Meta &meta=Meta());
      Data(const Vec<Scalar, 3>::Data &o, const std::string &n, Type id);
      ~Data();
      static Data *create(const std::string &name="", Type id = UNKNOWN, const Index numVertices = 0,
            const Meta &meta=Meta());

   V_DATA_END(Coords);
};

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Coords)

} // namespace vistle

#ifdef VISTLE_IMPL
#include "coords_impl.h"
#endif
#endif
