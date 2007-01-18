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
// $Id: DRawDataPoint.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

//template <Size D, typename Intensity = double, typename Position = float>
class DRawDataPoint1 	
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DRawDataPoint.h>
	using namespace OpenMS;
	typedef DRawDataPoint<1> DRawDataPoint1;
%End
	public:
	
	typedef DPosition1		PositionType;
	typedef double				IntensityType;
	DRawDataPoint1();
	DRawDataPoint1(const DRawDataPoint1& /p/);
	~DRawDataPoint1();
	DRawDataPoint1* clone() const;

	const IntensityType getIntensity() const;
	// IntensityType& getIntensity();
	void setIntensity(const IntensityType& /it/);

	const PositionType& getPosition() const;
	// PositionType& setPosition();
	void setPosition(const PositionType& /pos/);

	//DRawDataPoint& operator = (const DRawDataPoint& /rhs/)
	bool operator == (const DRawDataPoint1& /rhs/) const;
	bool operator != (const DRawDataPoint1& /rhs/) const;
	/*
	///
	class IntensityLess
	{
		public:
		IntensityLess(){}
		//IntensityLess(Size i){}
		~IntensityLess() {}
		
		bool operator () (const DRawDataPoint& a, const DRawDataPoint& b) const
	};
	
	///
	class PositionLess
	{
		public:
		PositionLess(Index i) : dimension_(i) {}
		PositionLess() : dimension_(-1) {}
		~PositionLess() {}

		bool operator () (const DRawDataPoint& a, const DRawDataPoint& b) const
	};
	*/
};
