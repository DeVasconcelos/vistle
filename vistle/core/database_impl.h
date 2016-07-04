#ifndef VISTLE_DATABASE_IMPL_H
#define VISTLE_DATABASE_IMPL_H

#include "scalars.h"
#include "archives.h"

#include <limits>

#include <boost/mpl/size.hpp>
#include "celltree_impl.h"

namespace vistle {

template<class Archive>
void DataBase::Data::serialize(Archive &ar, const unsigned int version) {
   ar & V_NAME("base_object", boost::serialization::base_object<Base::Data>(*this));
   ar & V_NAME("grid", grid);
}

} // namespace vistle

#endif
