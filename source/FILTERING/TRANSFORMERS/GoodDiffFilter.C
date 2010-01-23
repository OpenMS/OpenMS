// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  GoodDiffFilter::GoodDiffFilter()
    :FilterFunctor()
  {
		setName(GoodDiffFilter::getProductName());
    //values from kinter sherman


		// TODO from CHEMISTRY!
    aamass_.insert(make_pair(57.02 ,'G'));
    aamass_.insert(make_pair(71.04 ,'A'));
    aamass_.insert(make_pair(87.03 ,'S'));
    aamass_.insert(make_pair(97.05 ,'P'));
    aamass_.insert(make_pair(99.07 ,'V'));
    aamass_.insert(make_pair(101.05,'T'));
    aamass_.insert(make_pair(103.01,'C'));
    aamass_.insert(make_pair(113.08,'L')); //and also I, but the chars are fillers, anyway
    aamass_.insert(make_pair(114.04,'N'));
    aamass_.insert(make_pair(115.03,'D'));
    aamass_.insert(make_pair(128.06,'Q'));
    aamass_.insert(make_pair(128.09,'K'));
    aamass_.insert(make_pair(129.04,'E'));
    aamass_.insert(make_pair(131.04,'M'));
    aamass_.insert(make_pair(137.06,'H'));
    aamass_.insert(make_pair(147.07,'F'));
    aamass_.insert(make_pair(156.10,'R'));
    aamass_.insert(make_pair(163.06,'Y'));
    aamass_.insert(make_pair(186.06,'W'));
     
    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37, "Tolerance value as defined by Bern et al.");
		defaultsToParam_();
  }
  
  GoodDiffFilter::GoodDiffFilter( const GoodDiffFilter& source )
    :FilterFunctor(source),aamass_(source.aamass_)
  {
  }
  
  GoodDiffFilter& GoodDiffFilter::operator = (const GoodDiffFilter& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator=(source);
    	aamass_ = source.aamass_;
		}
    return *this;
  }

  GoodDiffFilter::~GoodDiffFilter()
  {
  }
}
