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
// $Id: ClusterNode.h,v 1.2 2006/03/29 13:06:18 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERNODE_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERNODE_H

#include <vector>
#include <list>

namespace OpenMS
{

  /**
  Node that represents the contents of a cluster of a MSMS-spectra<br>
  */
  class ClusterNode
  {
    friend class ClusterExperimentXMLHandler;
  public:
    typedef std::list<int>::iterator iterator;
    
    /** @brief consructor for merging two ClusterNodes <br> */
    ClusterNode(ClusterNode* left, ClusterNode* right); 
    
    /** @brief constructor for leaf <br>*/
    ClusterNode(int id,double parentmass = 0);
    
    /** @brief copy constructor <br> */
    ClusterNode(const ClusterNode& source);
    
    /** @brief standard constructor <br> */
    ClusterNode(); 
   
    /** @brief destructor <br> */
    ~ClusterNode();
    
    /** @brief assignment operator <br> */
    ClusterNode& operator=(const ClusterNode& source);
    
    /** @brief insert ClusterNode <br> */
    void insert(ClusterNode*);

    
    /**
    \return id if leaf, 0 else
    */
    int id() { return id_; }
    
    double getMaxParentMass() const {return max_parent_mass_;}
    double getMinParentMass() const {return min_parent_mass_;}
    
    /** @brief number of spectra in clster<br> */
    int size() {return childrenids_.size();}
    
    /** @brief readonly accessor for leafs */
    const std::list<int>& children() const;
  private:
    /** 
    ids of children
    */
    std::list<int> childrenids_;
    
    /**
    maximum parent mass of children
    */
    double max_parent_mass_;
    
    /**
    maximum parent mass of children
    */
    double min_parent_mass_;
    
    /** 
    > 0 if leaf
    */
    uint id_;
  };

}

#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTERNODE_H

