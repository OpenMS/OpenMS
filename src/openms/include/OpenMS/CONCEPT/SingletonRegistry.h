// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
