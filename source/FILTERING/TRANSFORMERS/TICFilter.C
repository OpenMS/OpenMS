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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <stdexcept>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cmath>
#include <ctime>

using namespace std;

namespace OpenMS
{
  TICFilter::TICFilter()
    :FilterFunctor()
  {
		setName(TICFilter::getProductName());
    defaults_.setValue("window", 5, "Windowing parameter which defines the windows size");
		defaultsToParam_();
  }

  TICFilter::TICFilter(const TICFilter& source)
    :FilterFunctor(source)
  {
  }

  TICFilter& TICFilter::operator = (const TICFilter& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator=(source);
		}
    return *this;
  }
  
  TICFilter::~TICFilter()
  {
  }

}
