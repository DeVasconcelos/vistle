#ifndef VISTLE_MAPDRAPE_MAPDRAPE_H
#define VISTLE_MAPDRAPE_MAPDRAPE_H

#include <vistle/module/module.h>
#include <vistle/core/vector.h>
#include <vistle/core/coords.h>
#ifdef MAPDRAPE
#include <vistle/module/resultcache.h>
#endif

#ifdef DISPLACE
#define MapDrape Displace
#endif

#ifdef DISPLACEVTKM
#define MapDrape DisplaceVtkm
#endif

using namespace vistle;

class MapDrape: public vistle::Module {
    static const unsigned NumPorts = 5;

public:
    MapDrape(const std::string &name, int moduleID, mpi::communicator comm);

private:
    bool compute() override;
    bool prepare() override;
    bool reduce(int timestep) override;

    Port *data_in[NumPorts], *data_out[NumPorts];

#ifdef MAPDRAPE
    StringParameter *p_mapping_from_, *p_mapping_to_;
    IntParameter *p_permutation;
    VectorParameter *p_offset;
    float offset[3];

    ResultCache<Coords::ptr> m_alreadyMapped;
#endif

#if defined(DISPLACE) || defined(DISPLACEVTKM)
    IntParameter *p_operation = nullptr;
    IntParameter *p_component = nullptr;
    FloatParameter *p_scale = nullptr;
#endif

#ifdef DISPLACEVTKM
    void applyDisplaceOperation(vistle::DataBase::const_ptr scalar, vistle::Coords::ptr coords, unsigned int c) const;
    void applyDisplaceWorklet(vistle::DataBase::const_ptr scalar, vistle::Coords::ptr coords) const;
#endif
};


#endif
