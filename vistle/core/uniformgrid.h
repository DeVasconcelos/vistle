//-------------------------------------------------------------------------
// UNIFORM GRID CLASS H
// *
// * Uniform Grid Container Object
//-------------------------------------------------------------------------
#ifndef UNIFORM_GRID_H
#define UNIFORM_GRID_H

#include "scalar.h"
#include "shm.h"
#include "structuredgridbase.h"
#include "export.h"

namespace vistle {

//-------------------------------------------------------------------------
// DECLARATION OF UNIFORMGRID
//-------------------------------------------------------------------------
class V_COREEXPORT UniformGrid : public StructuredGridBase {
   V_OBJECT(UniformGrid);

public:
   typedef StructuredGridBase Base;

   // constructor
   UniformGrid(Index xdim, Index ydim, Index zdim, const Meta &meta=Meta());

   // get/set functions for metadata
   Index getNumDivisions(int c) override { return (*d()->numDivisions)[c]; }
   Index getNumDivisions(int c) const override { return m_numDivisions[c]; }

   // get/set functions for shared memory members
   shm<Scalar>::array & min() { return *d()->min; }
   shm<Scalar>::array & max() { return *d()->max; }
   const Scalar * min() const { return m_min; }
   const Scalar * max() const { return m_max; }

private:
   // mutable pointers to ShmVectors
   mutable Index m_numDivisions[3];
   mutable Scalar m_min[3];
   mutable Scalar m_max[3];

   // data object
   V_DATA_BEGIN(UniformGrid);

   // each of the following variables represents a coordinate (by index, in order x, y, z)
   ShmVector<Index> numDivisions; //< number of divisions on each axis (1 more than number of cells)
   ShmVector<Scalar> min; //< coordinates of minimum grid point
   ShmVector<Scalar> max; //< coordinates of maximum grid point

   Data(const std::string & name, Index xdim, Index ydim, Index zdim, const Meta &meta=Meta());
   ~Data();
   static Data *create(Index xdim=0, Index ydim=0, Index zdim=0, const Meta &meta = Meta());

   V_DATA_END(UniformGrid);
};



BOOST_SERIALIZATION_ASSUME_ABSTRACT(UniformGrid)

} // namespace vistle

#ifdef VISTLE_IMPL
#include "uniformgrid_impl.h"
#endif

#endif /* UNIFORM_GRID_H */
