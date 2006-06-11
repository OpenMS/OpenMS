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
// $Id: NeutralLossMarker.h,v 1.3 2006/04/05 11:18:23 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

namespace OpenMS
{
  /**
  NeutralLossMarker marks peak pairs which could represent an ion an its neutral loss ( water, ammonia )<br>
  
  \param tolerance m/z tolerance
  */
  class NeutralLossMarker
    :public PeakMarker
  {
  public:
    /** @brief standard constructor <br> */
    NeutralLossMarker();

    /** @brief copy constructor <br> */
    NeutralLossMarker(const NeutralLossMarker& source);

    /** @brief destructor <br> */
    ~NeutralLossMarker();

    /** @brief assignment operator <br> */
    NeutralLossMarker& operator=(const NeutralLossMarker& source);

    static FactoryProduct* create() { return new NeutralLossMarker();}
    std::map<double,bool> operator()( MSSpectrum< DPeak<1> >&) const;
    String info() const;

		static const String getName()
		{
			return "NeutralLossMarker";
		}
  private:
    static const String info_;
  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H
