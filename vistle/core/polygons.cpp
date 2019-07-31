#include "polygons.h"

namespace vistle {

Polygons::Polygons(const Index numElements,
      const Index numCorners,
      const Index numVertices,
      const Meta &meta)
: Polygons::Base(Polygons::Data::create(numElements, numCorners, numVertices, meta))
{
    refreshImpl();
}

bool Polygons::isEmpty() {

   return Base::isEmpty();
}

bool Polygons::isEmpty() const {

   return Base::isEmpty();
}

void Polygons::refreshImpl() const {
}

bool Polygons::checkImpl() const {

   return true;
}

void Polygons::Data::initData() {
}

Polygons::Data::Data(const Polygons::Data &o, const std::string &n)
: Polygons::Base::Data(o, n)
{
   initData();
}

Polygons::Data::Data(const Index numElements, const Index numCorners,
                   const Index numVertices, const std::string & name,
                   const Meta &meta)
   : Polygons::Base::Data(numElements, numCorners, numVertices,
         Object::POLYGONS, name, meta)
{
   initData();
}


Polygons::Data * Polygons::Data::create(const Index numElements,
                            const Index numCorners,
                            const Index numVertices,
                            const Meta &meta) {

   const std::string name = Shm::the().createObjectId();
   Data *p = shm<Data>::construct(name)(numElements, numCorners, numVertices, name, meta);
   publish(p);

   return p;
}

V_OBJECT_TYPE(Polygons, Object::POLYGONS);
V_OBJECT_CTOR(Polygons);

} // namespace vistle
