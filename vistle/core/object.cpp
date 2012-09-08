#include <iostream>
#include <iomanip>

#include <limits.h>

#include "message.h"
#include "messagequeue.h"
#include "object.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

using namespace boost::interprocess;

namespace vistle {

template<int> size_t memorySize();

template<> size_t memorySize<4>() {

   return INT_MAX;
}

template<> size_t memorySize<8>() {

   return 68719476736;
}

Shm* Shm::s_singleton = NULL;

Shm::Shm(const int m, const int r, const size_t size,
         message::MessageQueue * mq)
   : m_moduleID(m), m_rank(r), m_objectID(0), m_messageQueue(mq) {

   m_shm = new managed_shared_memory(open_or_create, "vistle", size);
}

Shm::~Shm() {

   delete m_shm;
}

Shm & Shm::instance(const int moduleID, const int rank,
                    message::MessageQueue * mq) {

   if (!s_singleton)
      s_singleton = new Shm(moduleID, rank, memorySize<sizeof(void *)>(), mq);

   return *s_singleton;
}

Shm & Shm::instance() {

   if (!s_singleton)
      s_singleton = new Shm(-1, -1, memorySize<sizeof(void *)>(), NULL);

   return *s_singleton;
}

managed_shared_memory & Shm::getShm() {

   return *m_shm;
}

void Shm::publish(const shm_handle_t & handle) {

   if (m_messageQueue) {
      vistle::message::NewObject n(m_moduleID, m_rank, handle);
      m_messageQueue->getMessageQueue().send(&n, sizeof(n), 0);
   }
}

std::string Shm::createObjectID() {

   std::stringstream name;
   name << "Object_" << std::setw(8) << std::setfill('0') << m_moduleID
        << "_" << std::setw(8) << std::setfill('0') << m_rank
        << "_" << std::setw(8) << std::setfill('0') << m_objectID++;

   return name.str();
}

Object * Shm::getObjectFromHandle(const shm_handle_t & handle) {

   try {
      Object *object = static_cast<Object *>
         (m_shm->get_address_from_handle(handle));
      return object;

   } catch (interprocess_exception &ex) { }

   return NULL;
}

Object::Object(const Type type, const std::string & n,
               const int b, const int t): m_id(type), m_block(b), m_timestep(t) {

   size_t size = MIN(n.size(), 31);
   n.copy(m_name, size);
   m_name[size] = 0;
}

Object::~Object() {

}

Object::Type Object::getType() const {

   return m_id;
}

std::string Object::getName() const {

   return m_name;
}

int Object::getTimestep() const {

   return m_timestep;
}

int Object::getBlock() const {

   return m_block;
}

void Object::setTimestep(const int time) {

   m_timestep = time;
}

void Object::setBlock(const int blk) {

   m_block = blk;
}

Triangles::Triangles(const size_t numCorners, const size_t numVertices,
                     const std::string & name,
                     const int block, const int timestep)
   : Object(Object::TRIANGLES, name, block, timestep) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   m_x = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);
   m_y = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);
   m_z = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   m_cl = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numCorners, size_t(), alloc_inst_size_t);
}


Triangles * Triangles::create(const size_t numCorners,
                              const size_t numVertices,
                              const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   Triangles *t = static_cast<Triangles *>
      (Shm::instance().getShm().construct<Triangles>(name.c_str())[1](numCorners, numVertices, name, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(t);
   Shm::instance().publish(handle);
   */
   return t;
}

size_t Triangles::getNumCorners() const {

   return m_cl->size();
}

size_t Triangles::getNumVertices() const {

   return m_x->size();
}

Lines::Lines(const size_t numElements, const size_t numCorners,
             const size_t numVertices, const std::string & name,
             const int block, const int timestep)
   : Object(Object::LINES, name, block, timestep) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   m_x = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   m_y = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   m_z = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   m_el = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);

   m_cl = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numCorners, size_t(), alloc_inst_size_t);
}


Lines * Lines::create(const size_t numElements, const size_t numCorners,
                      const size_t numVertices,
                      const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   Lines *l = static_cast<Lines *>
      (Shm::instance().getShm().construct<Lines>(name.c_str())[1](numElements, numCorners, numVertices, name, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(l);
   Shm::instance().publish(handle);
   */
   return l;
}

size_t Lines::getNumElements() const {

   return m_el->size();
}

size_t Lines::getNumCorners() const {

   return m_cl->size();
}

size_t Lines::getNumVertices() const {

   return m_x->size();
}

Polygons::Polygons(const size_t numElements, const size_t numCorners,
                   const size_t numVertices, const std::string & name,
                   const int block, const int timestep)
   : Object(Object::POLYGONS, name, block, timestep) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   x = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   y = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   z = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   el = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);

   cl = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numCorners, size_t(), alloc_inst_size_t);
}


Polygons * Polygons::create(const size_t numElements,
                            const size_t numCorners,
                            const size_t numVertices,
                            const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   Polygons *p = static_cast<Polygons *>
      (Shm::instance().getShm().construct<Polygons>(name.c_str())[1](numElements, numCorners, numVertices, name, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(p);
   Shm::instance().publish(handle);
   */
   return p;
}

size_t Polygons::getNumElements() const {

   return el->size();
}

size_t Polygons::getNumCorners() const {

   return cl->size();
}

size_t Polygons::getNumVertices() const {

   return x->size();
}

UnstructuredGrid::UnstructuredGrid(const size_t numElements,
                                   const size_t numCorners,
                                   const size_t numVertices,
                                   const std::string & name,
                                   const int block, const int timestep)
   : Object(Object::UNSTRUCTUREDGRID, name, block, timestep) {

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   const allocator<size_t, managed_shared_memory::segment_manager>
      alloc_inst_size_t(Shm::instance().getShm().get_segment_manager());

   const allocator<char, managed_shared_memory::segment_manager>
      alloc_inst_char(Shm::instance().getShm().get_segment_manager());

   x = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   y = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   z = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numVertices, float(), alloc_inst_float);

   tl = Shm::instance().getShm().construct<boost::interprocess::vector<char, allocator<char, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, char(), alloc_inst_char);

   el = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, size_t(), alloc_inst_size_t);

   cl = Shm::instance().getShm().construct<boost::interprocess::vector<size_t, allocator<size_t, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numCorners, size_t(), alloc_inst_size_t);
}

UnstructuredGrid * UnstructuredGrid::create(const size_t numElements,
                                            const size_t numCorners,
                                            const size_t numVertices,
                                            const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   UnstructuredGrid *u = static_cast<UnstructuredGrid *>
      (Shm::instance().getShm().construct<UnstructuredGrid>(name.c_str())[1](numElements, numCorners, numVertices, name, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(u);
   Shm::instance().publish(handle);
   */
   return u;
}

size_t UnstructuredGrid::getNumElements() const {

   return el->size();
}

size_t UnstructuredGrid::getNumCorners() const {

   return cl->size();
}

size_t UnstructuredGrid::getNumVertices() const {

   return x->size();
}


Set::Set(const size_t numElements, const std::string & name,
         const int block, const int timestep)
   : Object(Object::SET, name, block, timestep) {

   const allocator<offset_ptr<Object>, managed_shared_memory::segment_manager>
      alloc_inst(Shm::instance().getShm().get_segment_manager());

   elements = Shm::instance().getShm().construct<boost::interprocess::vector<offset_ptr<Object>, allocator<offset_ptr<Object>, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(numElements, offset_ptr<Object>(), alloc_inst);
}


Set * Set::create(const size_t numElements,
                  const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   Set *p = static_cast<Set *>
      (Shm::instance().getShm().construct<Set>(name.c_str())[1](numElements, name, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(p);
   Shm::instance().publish(handle);
   */
   return p;
}

size_t Set::getNumElements() const {

   return elements->size();
}

Object * Set::getElement(const size_t index) const {

   if (index >= elements->size())
      return NULL;

   return (*elements)[index].get();
}

Geometry::Geometry(const std::string & name,
                   const int block, const int timestep)
   : Object(Object::GEOMETRY, name, block, timestep), geometry(NULL),
     colors(NULL), normals(NULL), texture(NULL) {

}


Geometry * Geometry::create(const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   Geometry *g = static_cast<Geometry *>
      (Shm::instance().getShm().construct<Geometry>(name.c_str())[1](name, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(g);
   Shm::instance().publish(handle);
   */
   return g;
}

Texture1D::Texture1D(const std::string & name, const size_t width,
                     const float mi, const float ma,
                     const int block, const int timestep)
   : Object(Object::TEXTURE1D, name, block, timestep), min(mi), max(ma) {

   const allocator<unsigned char, managed_shared_memory::segment_manager>
      alloc_inst(Shm::instance().getShm().get_segment_manager());

   pixels = Shm::instance().getShm().construct<boost::interprocess::vector<unsigned char, allocator<unsigned char, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(width * 4, char(), alloc_inst);

   const allocator<float, managed_shared_memory::segment_manager>
      alloc_inst_float(Shm::instance().getShm().get_segment_manager());

   coords = Shm::instance().getShm().construct<boost::interprocess::vector<float, allocator<float, managed_shared_memory::segment_manager> > >(Shm::instance().createObjectID().c_str())(alloc_inst_float);

}

Texture1D * Texture1D::create(const size_t width,
                              const float min, const float max,
                              const int block, const int timestep) {

   const std::string name = Shm::instance().createObjectID();
   Texture1D *tex = static_cast<Texture1D *>(Shm::instance().getShm().construct<Texture1D>(name.c_str())[1](name, width, min, max, block, timestep));

   /*
   shm_handle_t handle =
      Shm::instance().getShm().get_handle_from_address(tex);
   Shm::instance().publish(handle);
   */
   return tex;
}

size_t Texture1D::getNumElements() const {

   return coords->size();
}

size_t Texture1D::getWidth() const {

   return pixels->size() / 4;
}

template<> const Object::Type Vec<float>::s_type  = Object::VECFLOAT;
template<> const Object::Type Vec<int>::s_type    = Object::VECINT;
template<> const Object::Type Vec<char>::s_type   = Object::VECCHAR;
template<> const Object::Type Vec3<float>::s_type = Object::VEC3FLOAT;
template<> const Object::Type Vec3<int>::s_type   = Object::VEC3INT;
template<> const Object::Type Vec3<char>::s_type  = Object::VEC3CHAR;

} // namespace vistle
