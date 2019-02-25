#ifndef VISTLE_PARRENDMGR_H
#define VISTLE_PARRENDMGR_H

#include <IceT.h>
#include <IceTMPI.h>
#include "rhrcontroller.h"
#include "renderobject.h"

namespace vistle {

class RenderObject;
class Renderer;

void V_RENDEREREXPORT toIcet(IceTDouble *imat, const vistle::Matrix4 &vmat);

class V_RENDEREREXPORT ParallelRemoteRenderManager {
public:
   EIGEN_MAKE_ALIGNED_OPERATOR_NEW

   typedef void (*IceTDrawCallback)(const IceTDouble *proj, const IceTDouble *mv, const IceTFloat *bg, const IceTInt *viewport, IceTImage image);
   struct PerViewState;

   ParallelRemoteRenderManager(Renderer *module, IceTDrawCallback drawCallback);
   ~ParallelRemoteRenderManager();

   Port *outputPort() const;
   void connectionAdded(const Port *to);
   void connectionRemoved(const Port *to);
   Object::ptr getConfigObject();

   bool handleParam(const Parameter *p);
   bool prepareFrame(size_t numTimesteps);
   int timestep() const;
   size_t numViews() const;
   void setCurrentView(size_t i);
   void finishCurrentView(const IceTImage &img, int timestep);
   void finishCurrentView(const IceTImage &img, int timestep, bool lastView);
   bool finishFrame(int timestep);
   void getModelViewMat(size_t viewIdx, IceTDouble *mat) const;
   void getProjMat(size_t viewIdx, IceTDouble *mat) const;
   const PerViewState &viewData(size_t viewIdx) const;
   unsigned char *rgba(size_t viewIdx);
   float *depth(size_t viewIdx);
   void updateRect(size_t viewIdx, const IceTInt *viewport);
   void setModified();
   bool sceneChanged() const;
   bool isVariantVisible(const std::string &variant) const;
   void setLocalBounds(const Vector3 &min, const Vector3 &max);
   int rootRank() const {
      return m_displayRank==-1 ? 0 : m_displayRank;
   }
   void addObject(std::shared_ptr<RenderObject> ro);
   void removeObject(std::shared_ptr<RenderObject> ro);
   bool checkIceTError(const char *msg) const;

   Renderer *m_module;
   IceTDrawCallback m_drawCallback;
   int m_displayRank;
   RhrController m_rhrControl;
   IntParameter *m_continuousRendering;

   FloatParameter *m_delay;
   double m_delaySec;

   IntParameter *m_colorRank;
   Vector4 m_defaultColor;

   Vector3 localBoundMin, localBoundMax;

   size_t m_updateCount = -1;
   int m_updateBounds;
   int m_updateVariants;
   int m_updateScene;
   int m_doRender;
   size_t m_lightsUpdateCount;

   struct PerViewState {
       // synchronized across all ranks

      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      Matrix4 model;
      Matrix4 view;
      Matrix4 proj;
      std::vector<RhrServer::Light> lights;
      RhrServer::ViewParameters rhrParam;
      int width, height;

      PerViewState()
      : width(0)
      , height(0)
      {

         model.Identity();
         view.Identity();
         proj.Identity();
      }

      template<class Archive>
      void serialize(Archive &ar, const unsigned int version) {
         ar & width;
         ar & height;
         ar & model;
         ar & view;
         ar & proj;
         ar & lights;
      }
   };

   //! serializable description of one IceT tile - usable by boost::mpi
   struct DisplayTile {

      int rank;
      int x, y, width, height;

      DisplayTile()
      : rank(-1)
      , x(0)
      , y(0)
      , width(0)
      , height(0)
      {}

      template<class Archive>
      void serialize(Archive &ar, const unsigned int version) {
         ar & x;
         ar & y;
         ar & width;
         ar & height;
      }

   };

   //! state shared among all views
   struct GlobalState {
      int timestep = -1;
      int numTimesteps = 0;
      Vector3 bMin, bMax;

      GlobalState()
      : timestep(-1)
      , numTimesteps(0)
      {
      }

      template<class Archive>
      void serialize(Archive &ar, const unsigned int version) {
         ar & timestep;
         ar & numTimesteps;
      }
   };
   struct GlobalState m_state;

   std::vector<PerViewState, Eigen::aligned_allocator<PerViewState>> m_viewData; // synchronized from rank 0 to slaves
   std::vector<std::vector<unsigned char>> m_rgba;
   std::vector<std::vector<float>> m_depth;
   int m_currentView; //!< holds no. of view currently being rendered - not a problem as IceT is not reentrant anyway
   bool m_frameComplete; //!< track whether frame has been flushed to clients

   //! per view IceT state
   struct IceTData {

       bool ctxValid;
       int width, height; // dimensions of local tile
       IceTContext ctx;

       IceTData()
       : ctxValid(false)
       , width(0)
       , height(0)
       {
           ctx = 0;
       }
   };
   std::vector<IceTData> m_icet; // managed locally

   void updateVariants();
   RhrServer::VariantVisibilityMap m_clientVariants;
   Renderer::VariantMap m_localVariants;
};

}
#endif


