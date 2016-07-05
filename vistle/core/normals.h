#ifndef NORMALS_H
#define NORMALS_H


#include "vec.h"

namespace vistle {

class  V_COREEXPORT Normals: public Vec<Scalar, 3> {
   V_OBJECT(Normals);

   public:
   typedef Vec<Scalar,3> Base;

   Normals(const Index numNormals, const Meta &meta=Meta());

   Index getNumNormals() const;

   V_DATA_BEGIN(Normals);

      Data(const Index numNormals = 0,
            const std::string & name = "",
            const Meta &meta=Meta());
      Data(const Vec<Scalar, 3>::Data &o, const std::string &n);
      static Data *create(const Index numNormals = 0, const Meta &meta=Meta());

   V_DATA_END(Normals);
};

} // namespace vistle

V_OBJECT_DECLARE(vistle::Normals)

#ifdef VISTLE_IMPL
#include "normals_impl.h"
#endif

#endif
