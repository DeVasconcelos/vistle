﻿//
//This code is used for both IsoCut and IsoSurface!
//

#include <iostream>
#include <cmath>
#include <ctime>
#include <boost/mpi/collectives.hpp>
#include <core/message.h>
#include <core/object.h>
#include <core/unstr.h>
#include <core/vec.h>
#include "IsoSurface.h"
#include "Leveller.h"

MODULE_MAIN(IsoSurface)

using namespace vistle;

DEFINE_ENUM_WITH_STRING_CONVERSIONS(PointOrValue,
                                    (Point)
                                    (Value)
)

IsoSurface::IsoSurface(const std::string &shmname, const std::string &name, int moduleID)
   : Module(
#ifdef CUTTINGSURFACE
"cut through grids"
#else
"extract isosurfaces"
#endif
, shmname, name, moduleID) {

   setDefaultCacheMode(ObjectCache::CacheDeleteLate);
#ifdef CUTTINGSURFACE
   m_mapDataIn = createInputPort("data_in");
   addVectorParameter("point", "point on plane", ParamVector(0.0, 0.0, 0.0));
   addVectorParameter("vertex", "normal on plane", ParamVector(1.0, 0.0, 0.0));
   addFloatParameter("scalar", "distance to origin of ordinates", 0.0);
   m_option = addIntParameter("option", "option", Plane, Parameter::Choice);
   V_ENUM_SET_CHOICES(m_option, SurfaceOption);
   addVectorParameter("direction", "direction for variable Cylinder", ParamVector(0.0, 0.0, 0.0));
#else
   setReducePolicy(message::ReducePolicy::OverAll);
   m_isovalue = addFloatParameter("isovalue", "isovalue", 0.0);
   m_isopoint = addVectorParameter("isopoint", "isopoint", ParamVector(0.0, 0.0, 0.0));
   m_pointOrValue = addIntParameter("Interactor", "point or value interaction", Value, Parameter::Choice);
   V_ENUM_SET_CHOICES(m_pointOrValue, PointOrValue);

   createInputPort("data_in");
   m_mapDataIn = createInputPort("mapdata_in");
#endif    
   m_dataOut = createOutputPort("data_out");

   m_processortype = addIntParameter("processortype", "processortype", 0, Parameter::Choice);
   V_ENUM_SET_CHOICES(m_processortype, ThrustBackend);

   m_paraMin = m_paraMax = 0.f;
}

IsoSurface::~IsoSurface() {

}

bool IsoSurface::parameterChanged(const Parameter* param) {
#ifdef CUTTINGSURFACE
   switch(m_option->getValue()) {
      case Plane: {
         if(param->getName() == "point") {
            Vector vertex = getVectorParameter("vertex");
            Vector point = getVectorParameter("point");
            setFloatParameter("scalar",point.dot(vertex));
         }
         break;
      }
      case Sphere: {
         if(param->getName() == "point") {
            Vector vertex = getVectorParameter("vertex");
            Vector point = getVectorParameter("point");
            Vector diff = vertex - point;
            setFloatParameter("scalar",diff.norm());
         }break;
      }
      default:/*cylinders*/ {
         if(param->getName() == "option") {
            switch(m_option->getValue()) {

               case CylinderX:
                  setVectorParameter("direction",ParamVector(1,0,0));
                  break;
               case CylinderY:
                  setVectorParameter("direction",ParamVector(0,1,0));
                  break;
               case CylinderZ:
                  setVectorParameter("direction",ParamVector(0,0,1));
                  break;
            }
         }
         if(param->getName() == "point") {
            std::cerr<< "point" << std::endl;
            Vector vertex = getVectorParameter("vertex");
            Vector point = getVectorParameter("point");
            Vector direction = getVectorParameter("direction");
            Vector diff = direction.cross(vertex-point);
            setFloatParameter("scalar",diff.norm());
         }
         if(param->getName() == "scalar") {
            std::cerr<< "scalar" << std::endl;
            Vector vertex = getVectorParameter("vertex");
            Vector point = getVectorParameter("point");
            Vector direction = getVectorParameter("direction");
            Vector diff = direction.cross(vertex-point);
            float scalar = getFloatParameter("scalar");
            Vector newPoint = scalar*diff/diff.norm();
            setVectorParameter("point",ParamVector(newPoint[0], newPoint[1], newPoint[2]));
         }
         break;
      }
   }
#else
    if (param == m_isopoint) {
        setParameter(m_pointOrValue, (Integer)Point);
    } else if (param == m_isovalue) {
        setParameter(m_pointOrValue, (Integer)Value);
    }
#endif
   return true;
}

bool IsoSurface::prepare() {

   m_min = std::numeric_limits<Scalar>::max() ;
   m_max = -std::numeric_limits<Scalar>::max();
   return Module::prepare();
}

bool IsoSurface::reduce(int timestep) {

#ifndef CUTTINGSURFACE
   if (m_pointOrValue->getValue() == Point) {
       Scalar value = m_isovalue->getValue();
       int found = 0;
       Vector point = m_isopoint->getValue();
       for (size_t i=0; i<m_grids.size(); ++i) {
           if (m_datas[i]->getTimestep() == 0) {
               Index cell = m_grids[i]->findCell(point);
               if (cell != InvalidIndex) {
                   found = 1;
                   auto interpol = m_grids[i]->getInterpolator(cell, point);
                   value = interpol(m_datas[i]->x());
               }
           }
       }
       int numFound = boost::mpi::all_reduce(comm(), found, std::plus<int>());
       if (m_rank == 0 && numFound > 1) {
           sendWarning("found isopoint in %d blocks", numFound);
       }
       int valRank = found ? m_rank : m_size;
       valRank = boost::mpi::all_reduce(comm(), valRank, boost::mpi::minimum<int>());
       if (valRank < m_size) {
           boost::mpi::broadcast(comm(), value, valRank);
           setParameter(m_isovalue, (Float)value);
           setParameter(m_pointOrValue, (Integer)Point);
       }

       for (size_t i=0; i<m_grids.size(); ++i) {
           work(m_grids[i], m_datas[i], m_mapdatas[i]);
       }
       m_grids.clear();
       m_datas.clear();
       m_mapdatas.clear();
   }

   Scalar min, max;
   boost::mpi::all_reduce(comm(),
                          m_min, min, boost::mpi::minimum<Scalar>());
   boost::mpi::all_reduce(comm(),
                          m_max, max, boost::mpi::maximum<Scalar>());

   if (m_paraMin != (Float)min || m_paraMax != (Float)max)
      setParameterRange(m_isovalue, (Float)min, (Float)max);

   m_paraMax = max;
   m_paraMin = min;
#endif

   return Module::reduce(timestep);
}

bool IsoSurface::work(vistle::UnstructuredGrid::const_ptr gridS,
             vistle::Vec<vistle::Scalar>::const_ptr dataS,
             vistle::DataBase::const_ptr mapdata) {

   const int processorType = getIntParameter("processortype");
#ifdef CUTTINGSURFACE
   const Scalar isoValue = 0.0;
   int option = getIntParameter("option");
   const Vector vertex = getVectorParameter("vertex");
   const Vector point = getVectorParameter("point");
   const Vector direction = getVectorParameter("direction");
#else
   const Scalar isoValue = getFloatParameter("isovalue");
#endif

   Leveller l(gridS, isoValue, processorType
#ifdef CUTTINGSURFACE
              , option
#endif
              );

#ifdef CUTTINGSURFACE
   l.setCutData(vertex,point, direction);
#else
   l.setIsoData(dataS);
#endif
   if(mapdata){
      l.addMappedData(mapdata);
   };
   l.process();

#ifndef CUTTINGSURFACE
   auto minmax = dataS->getMinMax();
   if (minmax.first[0] < m_min)
      m_min = minmax.first[0];
   if (minmax.second[0] > m_max)
      m_max = minmax.second[0];
#endif

   Object::ptr result = l.result();
   DataBase::ptr mapresult = l.mapresult();
   if (result) {
#ifndef CUTTINGSURFACE
      result->copyAttributes(dataS);
#endif
      result->copyAttributes(gridS, false);
      if (mapdata && mapresult) {
         mapresult->copyAttributes(mapdata);
         mapresult->setGrid(result);
         addObject(m_dataOut, mapresult);
      }
#ifndef CUTTINGSURFACE
      else {
          addObject(m_dataOut, result);
      }
#endif
   }
   return true;
}

bool IsoSurface::compute() {

#ifdef CUTTINGSURFACE
   auto mapdata = expect<DataBase>(m_mapDataIn);
   if (!mapdata)
       return true;
   auto  gridS = UnstructuredGrid::as(mapdata->grid());
#else
   auto mapdata = accept<DataBase>(m_mapDataIn);
   auto dataS = expect<Vec<Scalar>>("data_in");
   if (!dataS)
      return true;
   if (dataS->guessMapping() != DataBase::Vertex) {
      sendError("need per-vertex mapping on data_in");
      return true;
   }
   auto  gridS = UnstructuredGrid::as(dataS->grid());
#endif
   if (!gridS) {
       sendError("grid required on input data");
       return true;
   }

#ifdef CUTTINGSURFACE
    return work(gridS, nullptr, mapdata);
#else
    if (m_pointOrValue->getValue() == Value) {
        return work(gridS, dataS, mapdata);
    } else {
        gridS->getCelltree();
        m_grids.push_back(gridS);
        m_datas.push_back(dataS);
        m_mapdatas.push_back(mapdata);
        return true;
    }
#endif
}
