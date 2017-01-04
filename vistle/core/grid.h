#ifndef VISTLE_GRID_H
#define VISTLE_GRID_H

#include "object.h"
#include "vector.h"
#include "database.h"
#include "geometry.h"
#include "export.h"

//#define INTERPOL_DEBUG

namespace vistle {

class V_COREEXPORT GridInterface: virtual public ElementInterface {
 public:

   enum FindCellFlags { NoFlags=0,
                        AcceptGhost=1,
                        ForceCelltree=2,
                        NoCelltree=4,
                      };

   virtual bool isGhostCell(Index elem) const = 0;
   virtual Index findCell(const Vector &point, int flags=NoFlags) const = 0;
   virtual bool inside(Index elem, const Vector &point) const = 0;
   virtual std::pair<Vector, Vector> cellBounds(Index elem) const = 0;
   virtual Scalar cellDiameter(Index elem) const = 0; //< approximate diameter of cell

   class Interpolator {
      std::vector<Scalar> weights;
      std::vector<Index> indices;

    public:
      Interpolator() {}
      Interpolator(std::vector<Scalar> &weights, std::vector<Index> &indices)
#if __cplusplus >= 201103L
      : weights(std::move(weights)), indices(std::move(indices))
#else
      : weights(weights), indices(indices)
#endif
      {
#ifndef NDEBUG
          check();
#endif
      }

      Scalar operator()(const Scalar *field) const {
         Scalar ret(0);
         for (size_t i=0; i<weights.size(); ++i)
            ret += field[indices[i]] * weights[i];
         return ret;
      }

      Vector3 operator()(const Scalar *f0, const Scalar *f1, const Scalar *f2) const {
         Vector3 ret(0, 0, 0);
         for (size_t i=0; i<weights.size(); ++i) {
            const Index ind(indices[i]);
            const Scalar w(weights[i]);
            ret += Vector3(f0[ind], f1[ind], f2[ind]) * w;
         }
         return ret;
      }

      bool check() const;
   };

   DEFINE_ENUM_WITH_STRING_CONVERSIONS(InterpolationMode, (First) // value of first vertex
         (Mean) // mean value of all vertices
         (Nearest) // value of nearest vertex
         (Linear) // barycentric/multilinear interpolation
         );

   virtual Interpolator getInterpolator(Index elem, const Vector &point, DataBase::Mapping mapping=DataBase::Vertex, InterpolationMode mode=Linear) const = 0;
   Interpolator getInterpolator(const Vector &point, DataBase::Mapping mapping=DataBase::Vertex, InterpolationMode mode=Linear) const {
       const Index elem = findCell(point);
       if (elem == InvalidIndex) {
           return Interpolator();
       }
       return getInterpolator(elem, point, mapping, mode);
   }
};

}
#endif
