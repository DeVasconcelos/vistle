#ifndef TRACER_BLOCKDATA_H
#define TRACER_BLOCKDATA_H

#include <vector>

#include <core/object.h>
#include <core/grid.h>
#include <core/lines.h>
#include <mutex>


class BlockData{

    friend class Particle;

private:
    vistle::Object::const_ptr m_grid;
    const vistle::GridInterface *m_gridInterface;
    vistle::Vec<vistle::Scalar, 3>::const_ptr m_vecfld;
    vistle::Vec<vistle::Scalar>::const_ptr m_scafld;
    vistle::DataBase::Mapping m_vecmap, m_scamap;
    const vistle::Scalar *m_vx, *m_vy, *m_vz, *m_p;

public:
    BlockData(vistle::Index i,
              vistle::Object::const_ptr grid,
              vistle::Vec<vistle::Scalar, 3>::const_ptr vdata,
              vistle::Vec<vistle::Scalar>::const_ptr pdata = nullptr);
    ~BlockData();

    const vistle::GridInterface *getGrid();
    vistle::Vec<vistle::Scalar, 3>::const_ptr getVecFld();
    vistle::DataBase::Mapping getVecMapping() const;
    vistle::Vec<vistle::Scalar>::const_ptr getScalFld();
    vistle::DataBase::Mapping getScalMapping() const;
    std::vector<vistle::Vec<vistle::Scalar, 3>::ptr> getIplVec();
    std::vector<vistle::Vec<vistle::Scalar>::ptr> getIplScal();
};
#endif
