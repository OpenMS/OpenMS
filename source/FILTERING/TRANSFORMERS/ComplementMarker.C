// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>

using namespace std;

namespace OpenMS
{
  ComplementMarker::ComplementMarker()
    : PeakMarker()
  {
		setName(ComplementMarker::getProductName());
    defaults_.setValue("tolerance", 1.0, "Tolerance value as defined by Bern et al.");
    defaults_.setValue("marks", 1, "How often a peak needs to be marked to be returned");
		defaultsToParam_();
  }

  ComplementMarker::ComplementMarker(const ComplementMarker& source)
    : PeakMarker(source)
  {
  }
  
  ComplementMarker::~ComplementMarker()
  {
  }

  ComplementMarker& ComplementMarker::operator = (const ComplementMarker& source)
  {
		if (this != &source)
		{
    	PeakMarker::operator=(source);
		}
    return *this;
  }

}
