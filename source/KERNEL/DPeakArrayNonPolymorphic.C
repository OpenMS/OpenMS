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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/KERNEL/DRawDataPoint.h>

namespace OpenMS
{	
	template <> 
  void DPeakArrayNonPolymorphic<1, DRawDataPoint<1> >::persistentWrite(PersistenceManager& /*pm*/, const char* /*name*/) const throw (Exception::Base)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
	template <> 
	void DPeakArrayNonPolymorphic<2, DRawDataPoint<2> >::persistentWrite(PersistenceManager& /*pm*/, const char* /*name*/) const throw (Exception::Base)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
	template <> 
	void DPeakArrayNonPolymorphic<3, DRawDataPoint<3> >::persistentWrite(PersistenceManager& /*pm*/, const char* /*name*/) const throw (Exception::Base)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}

	// Default instantiation of the class (forces complete compilation).
	DPeakArrayNonPolymorphic<1> default_dpeakarraynonpolymorphic_1;
	DPeakArrayNonPolymorphic<2> default_dpeakarraynonpolymorphic_2;
}
