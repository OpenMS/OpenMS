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


#include <OpenMS/DATASTRUCTURES/DataPoint.h>

namespace OpenMS
{
	
DataPoint::DataPoint()
{
	rt = 0;
	mz = 0;
	feature_id= 0;
	cluster_id = 0;
	charge=1;
	isotopes_per_peptide = 1;
	quality=0;
	cluster_size=0;
}
	
DataPoint::DataPoint(const DataPoint &copyin) : GridElement(copyin)
{
	feature_id = copyin.feature_id;
	intensities = copyin.intensities;
	cluster_id = copyin.cluster_id;
	charge = copyin.charge;
	isotopes_per_peptide = copyin.isotopes_per_peptide;
	quality = copyin.quality;
	mass_shifts = copyin.mass_shifts;
	cluster_size = copyin.cluster_size;
}

DataPoint& DataPoint::operator=(const DataPoint &rhs)
{
	this->rt = rhs.rt;
	this->mz = rhs.mz;
	this->intensities = rhs.intensities;
	this->feature_id = rhs.feature_id;
	this->cluster_id = rhs.cluster_id;
	this->mass_shifts = rhs.mass_shifts;
	this->cluster_size = rhs.cluster_size;
	return *this;
}

bool DataPoint::operator==(const DataPoint &rhs) const
{
	if( this->feature_id != rhs.feature_id) return false;
	if( this->rt != rhs.rt) return false;
	if( this->mz != rhs.mz) return false;
	if( this->intensities != rhs.intensities) return false;
	return true;
}


bool DataPoint::operator!=(const DataPoint &rhs) const
{
	return !(*this==rhs);
}
bool DataPoint::operator<(const DataPoint &rhs) const // required for built-in STL functions like sort
{
	if ( *this == rhs) return false;

	return this->feature_id < rhs.feature_id;
}

int DataPoint::getID()
{
	return feature_id;
}
}


