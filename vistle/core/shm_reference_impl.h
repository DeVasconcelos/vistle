#ifndef SHM_REFERENCE_IMPL_H
#define SHM_REFERENCE_IMPL_H

#include "archives_config.h"
#include "object.h"

namespace vistle {

template<class T>
template<class Archive>
void shm_ref<T>::save(Archive &ar) const {

    ar & V_NAME(ar, "shm_name", m_name);
    ar.template saveArray<typename T::value_type>(*this);
}

template<class T>
template<class Archive>
void shm_ref<T>::load(Archive &ar) {
   shm_name_t shmname;
   ar & V_NAME(ar, "shm_name", shmname);

   std::string arname = shmname.str();
   std::string name = ar.translateArrayName(arname);
   //std::cerr << "shm_ref: loading " << shmname << " for " << m_name << "/" << name << ", valid=" << valid() << std::endl;
   if (m_name.str().empty()) {
       unref();
       m_name = Shm::the().createArrayId(name);
   }
   if (name.empty() && !m_name.str().empty()) {
       ar.registerArrayNameTranslation(arname, m_name);
   }

   if (m_name.str() != name) {
       unref();
   }
   name = m_name.str();

   ObjectData *obj = ar.currentObject();
   if (obj) {
       //std::cerr << "obj " << obj->name << ": unresolved: " << name << std::endl;
       obj->unresolvedReference();
   }

   auto handler = ar.objectCompletionHandler();
   ar.template fetchArray<typename T::value_type>(arname, [this, arname, name, obj, handler]() -> void {
      auto ref = Shm::the().getArrayFromName<typename T::value_type>(name);
      //std::cerr << "shm_array: array completion handler: " << arname << " -> " << name << "/" << m_name << ", ref=" << ref << std::endl;
      if (!ref) {
          std::cerr << "shm_array: NOT COMPLETE: array completion handler: " << arname << " -> " << name << "/" << m_name << ", ref=" << ref << std::endl;
      }
      assert(ref);
      *this = ref;
      if (obj) {
          //std::cerr << "obj " << obj->name << ": RESOLVED: " << name << std::endl;
          obj->referenceResolved(handler);
      } else {
          std::cerr << "shm_array RESOLVED: " << name << ", but no handler" << std::endl;
      }
   });
   //std::cerr << "shm_array: first try: this=" << *this << ", ref=" << ref << std::endl;
}

template<class T>
T &shm_ref<T>::operator*() {
    return *m_p;
}

template<class T>
const T &shm_ref<T>::operator*() const {
    return *m_p;
}

template<class T>
T *shm_ref<T>::operator->() {
#ifdef NO_SHMEM
    return m_p;
#else
    return m_p.get();
#endif
}

template<class T>
const T *shm_ref<T>::operator->() const {
#ifdef NO_SHMEM
    return m_p;
#else
    return m_p.get();
#endif
}

template<class T>
int shm_ref<T>::refcount() const {
    if (m_p)
        return m_p->refcount();
    return -1;
}

template<class T>
const shm_name_t &shm_ref<T>::name() const {
    return m_name;
}

}
#endif
