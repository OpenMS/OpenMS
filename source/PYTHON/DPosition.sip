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
// $Id: DPosition.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//



//	template <Size D, typename CoordinateType = float>
class DPosition1
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DPosition.h>
	using namespace OpenMS;
	#include <strstream>
	typedef DPosition<1> DPosition1;
%End
	public:
	
	// !!!!
	typedef float CoordinateType;

	typedef CoordinateType ValueType;
	typedef CoordinateType* Iterator;
	typedef const CoordinateType* ConstIterator;
	typedef CoordinateType value_type;
	typedef CoordinateType& reference;
	typedef CoordinateType* pointer;
	typedef CoordinateType* iterator;
	typedef const CoordinateType* const_iterator;

	DPosition1();
	~DPosition1();

	DPosition1(const CoordinateType /x/);
	DPosition1(const DPosition1& /pos/);

	// ???? const CoordinateType& operator [] (Position /index/) const;
	// ???? CoordinateType& operator [] (Position /index/);

	CoordinateType operator * (const DPosition1& /point/) const throw();
	bool operator == (const DPosition1& /point/) const throw();
	bool operator != (const DPosition1& /point/) const throw();
	bool operator < (const DPosition1& /point/) const throw();
	bool operator > (const DPosition1& /point/) const throw();
	bool operator <= (const DPosition1& /point/) const throw();
	bool operator >= (const DPosition1& /point/) const throw();
	Size size() const;
	void clear();

	static const DPosition1 zero;
	static const DPosition1 min;
	static const DPosition1 max;

	/*
	ConstIterator begin() const throw();
	Iterator begin() throw();

	ConstIterator end() const throw();
	Iterator end() throw();
	*/


  SIP_PYOBJECT __str__();
%MethodCode
  std::strstream s;
  s << "DPosition1 { " << sipCpp->operator[](0) << " }\0";
  sipRes = PyString_FromString(s.str());
%End
  SIP_PYOBJECT __repr__();
%MethodCode
  std::strstream s;
  s << "DPosition1 { " << sipCpp->operator[](0) << " }\0";
  sipRes = PyString_FromString(s.str());
%End
  SIP_PYOBJECT __float__();
%MethodCode
  sipRes = PyFloat_FromDouble(sipCpp->operator [] (0));
%End
};

/* ????
template <Size D, typename CoordinateType>
const typename DPosition<D, CoordinateType>::DPosition DPosition<D, CoordinateType>::zero 
	= DPosition<D, CoordinateType>(0);

template <Size D, typename CoordinateType>
const typename DPosition<D, CoordinateType>::DPosition DPosition<D, CoordinateType>::min 
	= DPosition<D, CoordinateType>(std::numeric_limits<CoordinateType>::min());

template <Size D, typename CoordinateType>
const typename DPosition<D, CoordinateType>::DPosition DPosition<D, CoordinateType>::max 
	= DPosition<D, CoordinateType>(std::numeric_limits<CoordinateType>::max());
*/
