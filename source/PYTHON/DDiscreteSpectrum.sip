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
// $Id: DDiscreteSpectrum.sip,v 1.2 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// --------------------------------------------------------------------------
//


// template <Size D, typename Container = DPeakArray<D> > 
class DDiscreteSpectrum1
	: PersistenceManagerInterface, 
	 	MetaInfoInterface
{
%TypeHeaderCode
	#include <OpenMS/KERNEL/DDiscreteSpectrum.h>
	using namespace OpenMS;
	typedef DDiscreteSpectrum<1> DDiscreteSpectrum1;
%End
	public:

	// !!!!
	typedef DPeakArray1 Container;
	typedef Container ContainerType;
	typedef DPeak1 PeakType;
	typedef DSpectrumInfo1 InfoType;
	// typedef Container::iterator Iterator;
	// typedef Container::iterator iterator;
	// typedef Container::const_iterator ConstIterator;
	// typedef Container::const_iterator const_iterator;

	DDiscreteSpectrum1();
	DDiscreteSpectrum1(const DDiscreteSpectrum1& /rhs/);
	virtual ~DDiscreteSpectrum1();

	//	DDiscreteSpectrum& operator = (const DDiscreteSpectrum& rhs)

	Size size() const;

	const ContainerType& getContainer() const;
	// const ContainerType& getContainer();
	// ???? void setContainer(const ContainerType& /container/);

	// ???? const std::list<Identification>& identification() const {return identification_;}

	const InfoType& getInfo() const;
	// InfoType& getInfo();
	void setInfo(const InfoType& /info/);

	const DPosition1& getMin() const;
	const DPosition1& getMax() const;
	void updateMinMax();

	bool empty() const;

	/* ????	
	inline ConstIterator begin() const { return container_.begin(); }
	inline ConstIterator end() const { return container_.end(); }
	inline Iterator begin() { return container_.begin(); }
	inline Iterator end() { return container_.end(); }
	*/
};
