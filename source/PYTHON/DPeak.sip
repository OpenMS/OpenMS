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
// $Id: DPeak.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//


	//template <Size D,  
	//		typename Intensity = double, 
	//		typename Position = float, 
	//		typename Charge = SignedInt, 
	//		typename RValue = float, 
	//		typename Area = double, 
	//		typename FullWidthHalfMax = double>

class DPeak1 
	:	 DRawDataPoint1,
		 PersistenceManagerInterface, 
		 MetaInfoInterface
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DPeak.h>
	#include <strstream>
	using namespace OpenMS;
	typedef DPeak<1> DPeak1;
%End
	public:
	
	typedef DPosition1 PositionType;
	typedef double IntensityType;
	typedef DPosition1 WidthType;
	typedef double RValueType;
	typedef double AreaType;
	typedef double FullWidthHalfMaxType;
	typedef	SignedInt ChargeType;

	DPeak1(); 
	DPeak1(const DPeak1& /p/);
	virtual ~DPeak1();
	virtual DPeak1* clone() const;
	
	const RValueType getRValue() const;
	// RValueType& getRValue();
	void setRValue(RValueType& /value/);

	const AreaType getArea() const;
	// AreaType& getArea() const;
	void setArea(const AreaType& /area/);

	FullWidthHalfMaxType getFWHM() const;
	// FullWidthHalfMaxType& getFWHM() const;
	void setFWHM(const FullWidthHalfMaxType& /fwhm/);

	const DPosition1& getMaxPosition() const;
	// PositionType& getMaxPosition();
	void setMaxPosition(const DPosition1& /pos/);

	WidthType getWidth() const;
	// WidthType& getWidth();
	void setWidth(double /width/);

	ChargeType getCharge() const;
	// ChargeType& getCharge();
	void setCharge(Position /charge/);

	// DPeak1& operator = (const DPeak1& rhs);
	bool operator == (const DPeak1& /rhs/) const;
	// not necessary! virtual bool operator != (const DPeak1& /rhs/) const;

	// PersistenceManager interface
	void serialize(PersistenceManager& /f/);
	void deserialize_(const std::string& /attribute_name/, const std::string& /value/);

	//	Comparator classes.
/*
	class IntensityLess
	{
		public:
		 bool operator () (DPeak1 const& /left/, DPeak1 const& /right/) const;
		bool operator () (DPeak1 const * const /left/, DPeak1 const * const /right/ ) const;
	};
	
	/// 
	class PositionLess
	{
		public:

		PositionLess(Index /i/);
		PositionLess();
		~PositionLess();
		bool operator () (const DPeak1& /a/, const DPeak1& /b/);
	};
	
	///
	class WidthLess
	{
		public:
		WidthLess(Index i);
		WidthLess();
		~WidthLess();
		bool operator () (const DPeak1& /a/, const DPeak1& /b/);
	};	
*/

  // convert the DPeak to a string representation
	// ??? This should require additional code to wrap the
  // MetaInfo object of the peak as well, but then perhaps not.
  SIP_PYOBJECT __str__();
%MethodCode
	std::strstream s;
	s << "DPeak1 { " << sipCpp->getPosition()[0] << ", " << sipCpp->getIntensity()
		<< ", " << sipCpp->getWidth()[0] << ", " << sipCpp->getArea() << ", " << sipCpp->getFWHM()
		<< ", " << sipCpp->getRValue() << ", " << sipCpp->getCharge() << " }\0";
  sipRes = PyString_FromString(s.str());
%End
  SIP_PYOBJECT __repr__();
%MethodCode
	std::strstream s;
	s << "DPeak1 { " << sipCpp->getPosition()[0] << ", " << sipCpp->getIntensity()
		<< ", " << sipCpp->getWidth()[0] << ", " << sipCpp->getArea() << ", " << sipCpp->getFWHM()
		<< ", " << sipCpp->getRValue() << ", " << sipCpp->getCharge() << " }\0";
  sipRes = PyString_FromString(s.str());
%End
};
