#ifndef OBJECT_H
#define OBJECT_H


#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/string.hpp>

// include headers that implement an archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

#include "shm.h"

#include "serialize.h"

BOOST_CLASS_IMPLEMENTATION(vistle::shm<char>::string, boost::serialization::primitive_type)

namespace boost {
namespace serialization {

template<>
void access::destroy(const vistle::shm<char>::string *t);

template<>
void access::construct(vistle::shm<char>::string *t);

} // namespace serialization
} // namespace boost


namespace vistle {

namespace bi = boost::interprocess;

typedef bi::managed_shared_memory::handle_t shm_handle_t;
typedef char shm_name_t[32];

class Shm;

#define V_NAME(name, obj) \
   boost::serialization::make_nvp(name, (obj))

class Object {
   friend class Shm;
   friend class ObjectTypeRegistry;
#ifdef SHMDEBUG
   friend void Shm::markAsRemoved(const std::string &name);
#endif

public:
   typedef boost::shared_ptr<Object> ptr;
   typedef boost::shared_ptr<const Object> const_ptr;

   enum Type {
      UNKNOWN           = -1,
      VECFLOAT          =  0,
      VECINT            =  1,
      VECCHAR           =  2,
      VEC3FLOAT         =  3,
      VEC3INT           =  4,
      VEC3CHAR          =  5,
      TRIANGLES         =  6,
      LINES             =  7,
      POLYGONS          =  8,
      UNSTRUCTUREDGRID  =  9,
      SET               = 10,
      GEOMETRY          = 11,
      TEXTURE1D         = 12,
   };

   struct Info {
      uint64_t infosize;
      uint64_t itemsize;
      uint64_t offset;
      uint64_t type;
      uint64_t block;
      uint64_t timestep;
      virtual ~Info() {}; // for RTTI
   };

   virtual ~Object();

   shm_handle_t getHandle() const;

   Info *getInfo(Info *info = NULL) const;

   Type getType() const;
   std::string getName() const;

   int getBlock() const;
   int getTimestep() const;

   void setBlock(const int block);
   void setTimestep(const int timestep);

   void addAttribute(const std::string &key, const std::string &value = "");
   void setAttributeList(const std::string &key, const std::vector<std::string> &values);
   void copyAttributes(Object::const_ptr src, bool replace = true);
   bool hasAttribute(const std::string &key) const;
   std::string getAttribute(const std::string &key) const;
   std::vector<std::string> getAttributes(const std::string &key) const;

   void ref() const;
   void unref() const;
   int refcount() const;

   virtual void serialize(boost::archive::text_iarchive &ar) = 0;
   virtual void serialize(boost::archive::text_oarchive &ar) const = 0;
   virtual void serialize(boost::archive::binary_iarchive &ar) = 0;
   virtual void serialize(boost::archive::binary_oarchive &ar) const = 0;
   virtual void serialize(boost::archive::xml_iarchive &ar) = 0;
   virtual void serialize(boost::archive::xml_oarchive &ar) const = 0;

 protected:
   struct Data {
      Type type;
      shm_name_t name;
      int refcount;
      boost::interprocess::interprocess_mutex mutex;

      int block;
      int timestep;

      typedef shm<char>::string Attribute;
      typedef shm<char>::string Key;
      typedef shm<Attribute>::vector AttributeList;
      typedef std::pair<const Key, AttributeList> AttributeMapValueType;
      typedef shm<AttributeMapValueType>::allocator AttributeMapAllocator;
      typedef bi::map<Key, AttributeList, std::less<Key>, AttributeMapAllocator> AttributeMap;
      bi::offset_ptr<AttributeMap> attributes;
      void addAttribute(const std::string &key, const std::string &value = "");
      void setAttributeList(const std::string &key, const std::vector<std::string> &values);
      void copyAttributes(const Data *src, bool replace);
      bool hasAttribute(const std::string &key) const;
      std::string getAttribute(const std::string &key) const;
      std::vector<std::string> getAttributes(const std::string &key) const;

      Data(Type id, const std::string &name, int b, int t);
      ~Data();
      void ref();
      void unref();
      static Data *create(Type id, int b, int t);

      friend class boost::serialization::access;
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version);

      // not implemented
      Data &operator=(const Data &);
      Data();
      Data(const Data &);
   };

   Data *m_data;
 public:
   Data *d() const { return m_data; }
 protected:

   Object(Data *data);

   static void publish(const Data *d);
   static Object::ptr create(Data *);
 private:
   friend class boost::serialization::access;
   template<class Archive>
      void serialize(Archive &ar, const unsigned int version) {
         d()->serialize<Archive>(ar, version);
      }
};

} // namespace vistle

namespace boost {
namespace serialization {

template<>
void access::destroy(const vistle::Object::Data::AttributeList *t);

template<>
void access::construct(vistle::Object::Data::AttributeList *t);

template<>
void access::destroy(const vistle::Object::Data::AttributeMapValueType *t);

template<>
void access::construct(vistle::Object::Data::AttributeMapValueType *t);

} // namespace serialization
} // namespace boost


namespace vistle {

template<class Archive>
void Object::Data::serialize(Archive &ar, const unsigned int version) {
   ar & V_NAME("type", type);
   ar & V_NAME("block", block);
   ar & V_NAME("timestep", timestep);
   ar & V_NAME("attributes", *attributes);
}


class ObjectTypeRegistry {
   friend struct Object::Data;
   friend Object::ptr Object::create(Object::Data *);
   public:
   typedef Object::ptr (*CreateFunc)(Object::Data *d);
   typedef void (*DestroyFunc)(const std::string &name);

   template<class O>
   static void registerType(int id) {
      assert(s_typeMap.find(id) == s_typeMap.end());
      struct FunctionTable t = {
         O::createFromData,
         O::destroy,
      };
      s_typeMap[id] = t;
   }

   private:
   static CreateFunc getCreator(int id);
   static DestroyFunc getDestroyer(int id);

   struct FunctionTable {
      CreateFunc create;
      DestroyFunc destroy;
   };
   static const struct FunctionTable &getType(int id);
   typedef std::map<int, FunctionTable> TypeMap;
   static TypeMap s_typeMap;
};

//! declare a new Object type
#define V_OBJECT(Type) \
   public: \
   typedef boost::shared_ptr<Type> ptr; \
   typedef boost::shared_ptr<const Type> const_ptr; \
   typedef Type Class; \
   static boost::shared_ptr<const Type> as(boost::shared_ptr<const Object> ptr) { return boost::dynamic_pointer_cast<const Type>(ptr); } \
   static boost::shared_ptr<Type> as(boost::shared_ptr<Object> ptr) { return boost::dynamic_pointer_cast<Type>(ptr); } \
   static Object::ptr createFromData(Object::Data *data) { return Object::ptr(new Type(static_cast<Type::Data *>(data))); } \
   static void destroy(const std::string &name) { shm<Type::Data>::destroy(name); } \
   virtual void serialize(boost::archive::text_iarchive &tar) { \
      tar & *this; \
   } \
   virtual void serialize(boost::archive::text_oarchive &tar) const { \
      tar & *this; \
   } \
   virtual void serialize(boost::archive::binary_iarchive &bar) { \
      bar & *this; \
   } \
   virtual void serialize(boost::archive::binary_oarchive &bar) const { \
      bar & *this; \
   } \
   virtual void serialize(boost::archive::xml_iarchive &xar) { \
      xar & V_NAME("object", *this); \
   } \
   virtual void serialize(boost::archive::xml_oarchive &xar) const { \
      xar & V_NAME("object", *this); \
   } \
   protected: \
   struct Data; \
   Data *d() const { return static_cast<Data *>(m_data); } \
   Type(Data *data) : Base(data) {} \
   private: \
   friend class boost::serialization::access; \
   template<class Archive> \
      void serialize(Archive &ar, const unsigned int version) { \
         d()->template serialize<Archive>(ar, version); \
      } \
   friend boost::shared_ptr<const Object> Shm::getObjectFromHandle(const shm_handle_t &); \
   friend shm_handle_t Shm::getHandleFromObject(boost::shared_ptr<const Object>); \
   friend boost::shared_ptr<Object> Object::create(Object::Data *); \
   friend class ObjectTypeRegistry


//! register a new Object type (complex form, specify suffix for symbol names)
#define V_OBJECT_TYPE3(Type, suffix, id) \
   namespace { \
      class RegisterObjectType_##suffix { \
         public: \
                 RegisterObjectType_##suffix() { \
                    ObjectTypeRegistry::registerType<Type >(id); \
                 } \
      }; \
\
      static RegisterObjectType_##suffix registerObjectType_##suffix; \
   }

//! register a new Object type (simple form for non-templates, symbol suffix determined automatically)
#define V_OBJECT_TYPE(Type, id) \
   V_OBJECT_TYPE3(Type, Type, id)

} // namespace vistle

#endif
