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
// $Id: ClusterNode.C,v 1.1 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>


using namespace std;

namespace OpenMS 
{

  /**
  the leafs of the 2 trees are connected via the 
  (rightmost child of left).right = (leftmost child of right)<br>
  left and right are deleted <br>
  */
  ClusterNode::ClusterNode(ClusterNode* left, ClusterNode* right) 
    : childrenids_(),id_(0)
  {
    childrenids_.splice(childrenids_.end(),left->childrenids_);
    childrenids_.splice(childrenids_.end(),right->childrenids_);
    max_parent_mass_ = max(left->max_parent_mass_,right->max_parent_mass_);
    min_parent_mass_ = min(left->min_parent_mass_,right->min_parent_mass_);
    delete left;
    delete right;
  }

  ClusterNode::ClusterNode(const ClusterNode& source)
    :childrenids_(source.childrenids_),max_parent_mass_(source.max_parent_mass_),min_parent_mass_(source.min_parent_mass_), id_(source.id_)
  {
  }
 
  ClusterNode& ClusterNode::operator=(const ClusterNode& source)
  {
    childrenids_=source.childrenids_;
    max_parent_mass_ = source.max_parent_mass_;
    min_parent_mass_ = source.min_parent_mass_;
    id_=source.id_;
    return *this;
  }

  /** 
  newnode is deleted <br> 
  */
	void ClusterNode::insert(ClusterNode* newnode)
  {
    childrenids_.splice(childrenids_.end(),newnode->childrenids_);
    max_parent_mass_ = max(max_parent_mass_,newnode->max_parent_mass_);
    min_parent_mass_ = min(min_parent_mass_,newnode->min_parent_mass_);
    delete newnode;
	}
  
  //ids of all the children
  const list<int>& ClusterNode::children() const
  {
    return childrenids_;
  }

  ClusterNode::~ClusterNode()
  {
  }

  ClusterNode::ClusterNode(int id, double parentmass)
    : childrenids_(),max_parent_mass_(parentmass),min_parent_mass_(parentmass),id_(id)
  { 
    childrenids_.push_back(id);
  }

  ClusterNode::ClusterNode() 
    : childrenids_(),max_parent_mass_(0),min_parent_mass_(20000),id_(0)
  {
  }

}
