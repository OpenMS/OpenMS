// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: DSpectrumInfo.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//


// template <Size D>
class DSpectrumInfo1 
	: PersistenceManagerInterface
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DSpectrumInfo.h>
	using namespace OpenMS;
	typedef DSpectrumInfo<1> DSpectrumInfo1;
	typedef DPeak<1> DPeak1;
%End

	public:
	
	DSpectrumInfo1();
	typedef DPeak1 PeakType;

	const PeakType& getPrecursorPeak() const;
	// PeakType& getPrecursorPeak();
	void setPrecursorPeak(const PeakType& /peak/) const;

	double getNormalizedRetention() const;
	void setNormalizedRetention(double /rt/);
};
