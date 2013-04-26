/**************************************************************************\ 
 **                                                                        **
 **                                                                        **
 ** Description: Interface classes for application modules to the COVISE   **
 **              software environment                                      **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                                                                        **
 **                             (C)1997 RUS                                **
 **                Computer Center University of Stuttgart                 **
 **                            Allmandring 30                              **
 **                            70550 Stuttgart                             **
 ** Author:                                                                **
 ** Date:                                                                  **
\**************************************************************************/

#include "coRestraint.h"

#include <sstream>
#include <cstdio>

using namespace vistle;

//==========================================================================
//
//==========================================================================
coRestraint::coRestraint()
: changed(true)
, stringChanged(false)
{
}

//==========================================================================
//
//==========================================================================
coRestraint::~coRestraint()
{
}

//==========================================================================
//
//==========================================================================
void coRestraint::add(ssize_t mi, ssize_t ma)
{
   stringChanged = false;
   changed = true;
   min.push_back(mi);
   max.push_back(ma);
}


//==========================================================================
//
//==========================================================================
void coRestraint::add(ssize_t val)
{
   stringChanged = false;
   changed = true;
   min.push_back(val);
   max.push_back(val);
}


//==========================================================================
//
//==========================================================================
void coRestraint::add(const char *selection)
{
   stringChanged = false;
   changed = true;
   const char *c=selection;
   while(*c && (*c < '0' || *c >'9'))
      c++;
   while (*c)
   {
      int inc=0;
      ssize_t dumMax, dumMin;
      ssize_t numNumbers = sscanf(c,"%zd-%zd%n",&dumMin,&dumMax,&inc);
      if(numNumbers>0)
      {
         if(numNumbers==1)
         {
            dumMax=dumMin;
            if(inc == 0) // inc is 0 at least on windows if only one number is read
            {
               while(*c && (*c >= '0' && *c <='9'))
                  c++;
            }
         }
         min.push_back(dumMin);
         max.push_back(dumMax);
      }
      c += inc;
      while(*c && (*c < '0' || *c > '9'))
         c++;
   }
}


//==========================================================================
//
//==========================================================================
void coRestraint::clear()
{
   stringChanged = false;
   changed = true;
   min.clear();
   max.clear();
}


//==========================================================================
//
//==========================================================================
ssize_t coRestraint::lower() const
{
   ssize_t i=0, low;
   if (!min.empty())
      low = min[0];
   else
      return -1;
   while (i<min.size())
   {
      if ( (low>=min[i]) )
      {
         low = min[i];
      }
      i++;
   }
   return low;
}


//==========================================================================
//
//==========================================================================
ssize_t coRestraint::upper() const
{
   ssize_t i=0, up;
   if (!max.empty())
      up = max[0];
   else
      return -1;
   while (i<max.size())
   {
      if ( (up<=max[i]) )
      {
         up = max[i];
      }
      ++i;
   }
   return up;
}


//==========================================================================
//
//==========================================================================
ssize_t coRestraint::operator ()(ssize_t val) const
{
   ssize_t i=0;
   while (i<min.size())
   {
      if ( (val>=min[i]) && (val<=max[i]) )
         return 1;
      i++;
   }
   return 0;
}


//==========================================================================
//
//==========================================================================
ssize_t coRestraint::get(ssize_t val, ssize_t &group) const
{
   group=0;
   while (group<min.size())
   {
      if ( (val>=min[group]) && (val<=max[group]) )
         return 1;
      group++;
   }
   return 0;
}

//==========================================================================
//
//==========================================================================

const std::string &coRestraint::getRestraintString() const
{
   if (!stringChanged)
   {
      stringChanged = true;
      restraintString = getRestraintString(getValues());
   }
   return restraintString;
}

const std::string coRestraint::getRestraintString(std::vector<ssize_t> sortedValues) const
{
   std::ostringstream restraintStream;
   ssize_t old=-1, size=sortedValues.size();
   bool started=false, firsttime=true;
   if (size == 0)
      return "";
   for(ssize_t i=0; i<size; ++i)
   {
      ssize_t actual = sortedValues[i];
      if (firsttime)
      {
         firsttime = false;
         restraintStream << actual;
         old=actual;
         continue;
      }
      else if ( actual == old+1 && i<size-1)
      {
         if (!started)
         {
            restraintStream << "-";
            started = true;
         }
         old = actual;
         continue;
      }
      else if (started)
      {
         started = false;
         restraintStream << old << ", " << actual;
      }
      else 
      {
         restraintStream << ", " << actual;
      }
      old = actual;
   }
   return restraintStream.str();
}

//==========================================================================
//
//==========================================================================
// function returns vector containing all integer indices
// that are specified by the string added to this coRestraint object
//
// returns an empty vector, if the evaluation of char array is 
// not successful
// 
const std::vector<ssize_t> &coRestraint::getValues() const
{
   if (changed)
   {
      changed = false;
      values.clear();
      //getting the indices
      ssize_t counter = lower();
      ssize_t limit = upper();
      if (limit == -1 || counter == -1)
      {
         values.clear();
      }
      else
      {
         while( counter <= limit )
         {
            if( operator()(counter) )
            {
               values.push_back(counter);
            }
            ++counter;
         }
      }
   }

   return values;
}

