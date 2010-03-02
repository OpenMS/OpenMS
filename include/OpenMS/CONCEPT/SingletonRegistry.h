	// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl, Chris Bielow  $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_SINGLETONREGISTRY_H
#define OPENMS_CONCEPT_SINGLETONREGISTRY_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <map>
#include <iostream>

namespace OpenMS
{
  class String;
	class FactoryBase;

  /**
  	@brief Holds pointers to unique instance of a singleton factory.
 		
		@note: NEVER(!) include this file anywhere (except for the SingletonRegistry.C)! :D
		
 		@ingroup Concept
  */
  class OPENMS_DLLAPI SingletonRegistry
  {
    friend class singletonsNeedNoFriends; //some versions of gcc would warn otherwise

  private:
    /// Function signature of creator function 
    typedef std::map<String, FactoryBase*> Map;
    typedef Map::const_iterator MapIterator;

    /// destructor 
    virtual ~SingletonRegistry(){}

    /// C'Tor
    SingletonRegistry(){}

    /// singleton access to SingletonRegistry 
    static SingletonRegistry* instance_()
    {
      if (!singletonRegistryInstance_)
			{
				singletonRegistryInstance_ = new SingletonRegistry();
      }
      return singletonRegistryInstance_;
    }

  public:
    
		/// return DefaultParamHandler according to unique identifier @p name  
    static FactoryBase* getFactory(const String& name)
    {
    	MapIterator it = instance_()->inventory_.find(name);
      if (it != instance_()->inventory_.end())
			{
				return it->second;
			}
      else 
			{
      	throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This Factory is not registered with SingletonRegistry!",name.c_str());
			}
    }
    
    /**
    	@brief register new concrete Factory 
     
       \param name unique name for Factory of certain type
       \param instance pointer to this Factory
    */
    static void registerFactory(const String& name, FactoryBase* instance)
    {
      instance_()->inventory_[name] = instance;
    }
		
		/// Returns if a factory is registered
		static bool isRegistered(String name)
		{
      if (instance_()->inventory_.find(name) != instance_()->inventory_.end())
			{
				return true;
			}
			return false;
		}
		
  private:

    Map inventory_;
    static SingletonRegistry* singletonRegistryInstance_;
  };
  
  
}
#endif //OPENMS_CONCEPT_SINGLETONREGISTRY_H
