// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_FACTORY_H
#define OPENMS_CONCEPT_FACTORY_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/FactoryBase.h>
#include <OpenMS/CONCEPT/SingletonRegistry.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <iostream>
#include <typeinfo>

namespace OpenMS
{
  class String;

  /**
  	@brief Returns FactoryProduct* based on the name of the desired concrete FactoryProduct
 		
 		Every factory product base class T has to implement the static function registerChildren that registers all classes S derived from T at Factory<T>.

		Every class S derived from T has to implement the function "static T* create()" which is going to be registered at Factory<T>.
		Additionally the function "static String getProductName()" is required, which returns the name the class is registered by.
 		
 		@ingroup Concept
  */
  template <typename FactoryProduct>
  class Factory
	: public FactoryBase
  {
    friend class singletonsNeedNoFriends; //some versions of gcc would warn otherwise

  private:
    /// Function signature of creator function 
    typedef FactoryProduct* (*FunctionType)();
    typedef std::map<String, FunctionType> Map;
    typedef typename Map::const_iterator MapIterator;
		typedef Factory<FactoryProduct> FactoryType;
		
    /// Destructor 
    virtual ~Factory(){}

    /// Constructor
    Factory()
    {
    }

    /// singleton access to Factory 
    static Factory* instance_()
    {
			if (!instance_ptr_)
			{
				// name of this Factory
				String myName = typeid(FactoryType).name();
				
				//check if an instance of this kind of Factory already registered
				if (!SingletonRegistry::isRegistered(myName))
				{
					// if not registered yet ... add it
					instance_ptr_ = new Factory();
					// now, attention as ORDER of commands is important here:
					// first register the Factory
					SingletonRegistry::registerFactory(myName, instance_ptr_);
					// because this call, might use another instance of this factory, but we want the other instance to register the children with "US"
					FactoryProduct::registerChildren();
				}
				else
				{
					// get instance of this factory from registry
					instance_ptr_ = static_cast<FactoryType*> (SingletonRegistry::getFactory(myName));			
				}
			}
      return instance_ptr_;
    }

  public:
		
    /// return FactoryProduct according to unique identifier @p name  
    static FactoryProduct* create(const String& name)
    {
    	MapIterator it = instance_()->inventory_.find(name);
      if (it != instance_()->inventory_.end())
			{
				return (*(it->second))();
			}
      else 
			{
      	throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "This FactoryProduct is not registered!",name.c_str());
			}
    }
    
    /**
    	@brief register new concrete FactoryProduct 
     
       @param name unique name for concrete FactoryProduct
       @param creator default constructor for concrete FactoryProduct 
    */
    static void registerProduct(const String& name, const FunctionType creator)
    {
      instance_()->inventory_[name] = creator;
    }
		
		/// Returns if a factory product is registered
		static bool isRegistered(const String& name)
		{
      if (instance_()->inventory_.find(name) != instance_()->inventory_.end())
			{
				return true;
			}
			return false;
		}

		/// Returns a list of registered products
		static std::vector<String> registeredProducts()
		{
			std::vector<String> list;
      for (MapIterator it = instance_()->inventory_.begin(); it!=instance_()->inventory_.end(); ++it)
      {
				list.push_back(it->first);
			}
			return list;
		}
		
  private:

    Map inventory_;
    static Factory* instance_ptr_;
  };
  
  template<typename FactoryProduct> Factory<FactoryProduct>* Factory<FactoryProduct>::instance_ptr_ = 0;
    
}
#endif //OPENMS_CONCEPT_FACTORY_H
