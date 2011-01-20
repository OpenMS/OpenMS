// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>
#include <OpenMS/CONCEPT/Factory.h>

using namespace std;

namespace OpenMS
{
  ClusterFunctor::ClusterFunctor()
	{
	}
  
	ClusterFunctor::ClusterFunctor(const ClusterFunctor& /*source*/)
	{
	}
 
	ClusterFunctor::~ClusterFunctor()
	{
	}
 
	ClusterFunctor& ClusterFunctor::operator=(const ClusterFunctor& /*source*/)
	{
		return *this;
	}
 
	void ClusterFunctor::registerChildren()
	{
		Factory<ClusterFunctor>::registerProduct(SingleLinkage::getProductName(), &SingleLinkage::create);
		Factory<ClusterFunctor>::registerProduct(CompleteLinkage::getProductName(), &CompleteLinkage::create);
		Factory<ClusterFunctor>::registerProduct(AverageLinkage::getProductName(), &AverageLinkage::create);
	}
	
	ClusterFunctor::InsufficientInput::InsufficientInput(const char* file, int line, const char* function, const char* message) throw()
          : BaseException(file, line, function, "ClusterFunctor::InsufficentInput", message)
	{
	} 

	ClusterFunctor::InsufficientInput::~InsufficientInput() throw()
	{
	}

}
