// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow  $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <map>

namespace OpenMS
{
  class String;
  class FactoryBase;

  /**
    @brief Holds pointers to unique instance of a singleton factory.

        @note: NEVER(!) include this file anywhere (except for the SingletonRegistry.cpp)! :D

        @ingroup Concept
  */
  class OPENMS_DLLAPI SingletonRegistry
  {
    friend class singletonsNeedNoFriends; //some versions of gcc would warn otherwise

private:
    /// Function signature of creator function
    typedef std::map<String, FactoryBase *> Map;
    typedef Map::const_iterator MapIterator;

    /// destructor
    virtual ~SingletonRegistry(){}

    /// C'Tor
    SingletonRegistry(){}

    /// singleton access to SingletonRegistry
    static SingletonRegistry * instance_()
    {
      if (!singletonRegistryInstance_)
      {
        singletonRegistryInstance_ = new SingletonRegistry();
      }
      return singletonRegistryInstance_;
    }

public:

    /// return DefaultParamHandler according to unique identifier @p name
    static FactoryBase * getFactory(const String & name)
    {
      MapIterator it = instance_()->inventory_.find(name);
      if (it != instance_()->inventory_.end())
      {
        return it->second;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "This Factory is not registered with SingletonRegistry!", name.c_str());
      }
    }

    /**
        @brief register new concrete Factory

       \param name unique name for Factory of certain type
       \param instance pointer to this Factory
    */
    static void registerFactory(const String & name, FactoryBase * instance)
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
    static SingletonRegistry * singletonRegistryInstance_;
  };


}
