#include "lines.h"
#include "lines_impl.h"
#include "archives.h"

namespace vistle {

Lines::Lines(const Index numElements, const Index numCorners,
                      const Index numVertices,
                      const Meta &meta)
   : Lines::Base(Lines::Data::create("", numElements, numCorners,
            numVertices, meta))
{
    refreshImpl();
}

void Lines::refreshImpl() const {
}

bool Lines::isEmpty() {

   return Base::isEmpty();
}

bool Lines::isEmpty() const {

   return Base::isEmpty();
}

bool Lines::checkImpl() const {

   return true;
}

void Lines::Data::initData() {
}

Lines::Data::Data(const Data &other, const std::string &name)
: Lines::Base::Data(other, name)
{
   initData();
}

Lines::Data::Data(const Index numElements, const Index numCorners,
             const Index numVertices, const std::string & name,
             const Meta &meta)
   : Lines::Base::Data(numElements, numCorners, numVertices,
         Object::LINES, name, meta)
{
   initData();
}


Lines::Data * Lines::Data::create(const std::string &objId, const Index numElements, const Index numCorners,
                      const Index numVertices,
                      const Meta &meta) {

   const std::string name = Shm::the().createObjectId(objId);
   Data *l = shm<Data>::construct(name)(numElements, numCorners, numVertices, name, meta);
   publish(l);

   return l;
}

V_OBJECT_TYPE(Lines, Object::LINES)
V_OBJECT_CTOR(Lines)
V_OBJECT_IMPL(Lines)

} // namespace vistle
