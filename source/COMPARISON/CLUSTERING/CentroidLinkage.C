// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/COMPARISON/CLUSTERING/CentroidLinkage.h>

namespace OpenMS
{
CentroidLinkage::CentroidLinkage(DoubleReal rt_scaling_) : ClusteringMethod(rt_scaling_)
{
}
DoubleReal CentroidLinkage::getDistance(DataSubset& subset1,DataSubset& subset2)
{
	return sqrt(pow(subset1.rt-subset2.rt,2)*rt_scaling*rt_scaling+pow(subset1.mz-subset2.mz,2));
}
DoubleReal CentroidLinkage::getDistance(DataPoint& point1,DataPoint& point2)
{
	return sqrt(pow(point1.rt-point2.rt,2)*rt_scaling*rt_scaling+pow(point1.mz-point2.mz,2));
}
}

