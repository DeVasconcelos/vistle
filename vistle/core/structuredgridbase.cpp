//-------------------------------------------------------------------------
// STRUCTURED GRID OBJECT BASE CLASS CPP
// *
// * Base class for Structured Grid Objects
//-------------------------------------------------------------------------

#include "structuredgridbase.h"

namespace vistle {

#ifndef _WIN32
constexpr std::array<Index,8> StructuredGridBase::HexahedronIndices[3];
#endif


// IS GHOST CELL CHECK
//-------------------------------------------------------------------------
bool StructuredGridBase::isGhostCell(Index elem) const {
      if (elem == InvalidIndex)
          return false;

      Index dims[3];
      std::array<Index,3> cellCoords;

      for (int c=0; c<3; ++c) {
          dims[c] = getNumDivisions(c);
      }

      cellCoords = cellCoordinates(elem, dims);

      for (int c=0; c<3; ++c) {
          if (cellCoords[c] < getNumGhostLayers(c, Bottom)
                  || cellCoords[c]+getNumGhostLayers(c, Top)+1 >= dims[c]) {
              return true;
          }
      }

      return false;

}

Scalar StructuredGridBase::cellDiameter(Index elem) const {

    auto bounds = cellBounds(elem);
    return (bounds.second-bounds.first).norm();
}

Vector StructuredGridBase::cellCenter(Index elem) const
{
    Index dims[3];
    for (int c=0; c<3; ++c) {
        dims[c] = getNumDivisions(c);
    }

    auto verts = cellVertices(elem, dims);
    Vector center(0);
    for (auto v: verts) {
        center += getVertex(v);
    }
    center *= 1./8;
    return center;

}

std::vector<Index> StructuredGridBase::getNeighborElements(Index elem) const {

    std::vector<Index> elems;
    if (elem == InvalidIndex)
        return elems;

    const Index dims[3] = { getNumDivisions(0), getNumDivisions(1), getNumDivisions(2) };
    const auto coords = cellCoordinates(elem, dims);
    for (int d=0; d<3; ++d) {
        auto c = coords;
        if (coords[d] >= 1) {
            c[d] = coords[d]-1;
            elems.push_back(cellIndex(c[0], c[1], c[2], dims));
        }
        if (coords[d] < dims[d]-2) {
            c[d] = coords[d]+1;
            elems.push_back(cellIndex(c[0], c[1], c[2], dims));
        }
    }

    return elems;
}

} // namespace vistle
