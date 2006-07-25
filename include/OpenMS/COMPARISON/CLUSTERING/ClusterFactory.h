// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERFACTORY_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERFACTORY_H

#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <map>
#include <vector>

namespace OpenMS
{

  /**
  returns FactoryProduct* based on the name of a FactoryProduct<br>
  New FactoryProducts should be registered in init(), but can also be registered temporarily with registerfp()
  
  @todo replace by Factory
  */
  class ClusterFactory
  {
    friend class singletonsNeedNoFriends; //some versions of gcc would warn otherwise
  public:
    /** @brief return FactoryProduct with unique identifier <i>name</i> <br> */
   	FactoryProduct* create(String name) const;

    /** @brief temporary registration <br> */
    void registerfp(String name, FactoryProduct*(*vp)() );

    /** @brief singleton access to ClusterFactory <br> */
    static ClusterFactory* instance();

    /** @brief permanent registration <br> */
    void init();

    /** @brief use instead of destructor */
    void destroy() { delete instancep_;instancep_ = 0; }

    /** @brief return registered FactoryProducts <br> */
    std::vector<String> catalogue(String type = "FactoryProduct") const;

    /** @brief copy a FactoryProduct <br> */
    FactoryProduct* duplicate(const FactoryProduct* tmplate ) const;

  private:
    /** @brief destructor <br> */
    virtual ~ClusterFactory();

    /** @brief create with instance <br> */
    ClusterFactory();

    /** @brief create with instance <br> */
    ClusterFactory(const ClusterFactory&);

    /** @brief create with instance <br> */
    ClusterFactory& operator=(const ClusterFactory& source);

    std::map<String, FactoryProduct*(*)()> inventory_;

    static ClusterFactory* instancep_ ;
  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTERFACTORY_H
