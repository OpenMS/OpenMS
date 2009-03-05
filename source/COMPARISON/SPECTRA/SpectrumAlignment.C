// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

using namespace std;

namespace OpenMS
{
  SpectrumAlignment::SpectrumAlignment()
		: DefaultParamHandler("SpectrumAlignment")
  {
		defaults_.setValue("tolerance", 0.3, "Defines the absolut (in Da) or relative (in ppm) tolerance");
		defaults_.setValue("is_relative_tolerance", "false", "If true, the 'tolerance' is interpreted as ppm-value");
		defaults_.setValidStrings("is_relative_tolerance", StringList::create("true,false"));
		defaultsToParam_();
  }

  SpectrumAlignment::SpectrumAlignment(const SpectrumAlignment& source)
		: DefaultParamHandler(source)
  {
  }

  SpectrumAlignment::~SpectrumAlignment()
  {
  }

  SpectrumAlignment& SpectrumAlignment::operator = (const SpectrumAlignment& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator = (source);
		}
    return *this;
  }
	
}
