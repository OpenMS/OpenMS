// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
  class Factory :
    public FactoryBase
  {
    friend class singletonsNeedNoFriends; //some versions of gcc would warn otherwise

private:
    /// Function signature of creator function
    typedef FactoryProduct * (*FunctionType)();
    typedef std::map<String, FunctionType> Map;
    typedef typename Map::const_iterator MapIterator;
    typedef Factory<FactoryProduct> FactoryType;

    /// Destructor
    ~Factory() override{}

    /// Constructor
    Factory()
    {
    }

    /// singleton access to Factory
    static Factory * instance_()
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
          instance_ptr_ = static_cast<FactoryType *>(SingletonRegistry::getFactory(myName));
        }
      }
      return instance_ptr_;
    }

public:

    /// return FactoryProduct according to unique identifier @p name
    static FactoryProduct * create(const String & name)
    {
      MapIterator it = instance_()->inventory_.find(name);
      if (it != instance_()->inventory_.end())
      {
        return (*(it->second))();
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "This FactoryProduct is not registered!", name.c_str());
      }
    }

    /**
        @brief register new concrete FactoryProduct

       @param name unique name for concrete FactoryProduct
       @param creator default constructor for concrete FactoryProduct
    */
    static void registerProduct(const String & name, const FunctionType creator)
    {
      instance_()->inventory_[name] = creator;
    }

    /// Returns if a factory product is registered
    static bool isRegistered(const String & name)
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
      for (MapIterator it = instance_()->inventory_.begin(); it != instance_()->inventory_.end(); ++it)
      {
        list.push_back(it->first);
      }
      return list;
    }

private:

    Map inventory_;
    static Factory * instance_ptr_;
  };

  template <typename FactoryProduct>
  Factory<FactoryProduct> * Factory<FactoryProduct>::instance_ptr_ = nullptr;

}
#endif //OPENMS_CONCEPT_FACTORY_H
