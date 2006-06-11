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
// $Id: PeakMarker.h,v 1.6 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H

#include <vector>
#include <map>
#include <string>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{
  /**
  a PeakMarker marks peaks that seem to fulfil some criterion<br>
  */
  class PeakMarker
    : public FactoryProduct
  {
  public:
    /** @brief standard constructor <br> */
    PeakMarker() ;

    /** @brief copy constructor <br> */
    PeakMarker(const PeakMarker& source);

    /** @brief destructor <br> */
    virtual ~PeakMarker() {}

    /** @brief assignment operator <br> */
    PeakMarker& operator=(const PeakMarker& source);
    virtual std::map<double,bool> operator()( MSSpectrum< DPeak<1> >&) const = 0;
  };

}
#endif // OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H
