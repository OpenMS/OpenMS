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
// $Id: ThresholdMower.h,v 1.4 2006/04/05 11:18:23 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/MowerFunctor.h>

namespace OpenMS
{
  /**
  ThresholdMover removes all peaks below a Threshold<br>
  
  \param threshold: the threshold
  for comparable results we suggest normalizing (for example with Normalizer) all
  Spectra first
  */
  class ThresholdMower
    :public MowerFunctor
  {
  public:
    /** @brief standard constructor <br> */
    ThresholdMower();

    /** @brief copy constructor <br> */
    ThresholdMower(const ThresholdMower& source);

    /** @brief destructor <br> */
    ~ThresholdMower();

    /** @brief assignment operator <br> */
    ThresholdMower& operator=(const ThresholdMower& source);

    static FactoryProduct* create() { return new ThresholdMower();}
    void operator()(MSSpectrum< DPeak<1> >&) const;
    String info() const;

		static const String getName()
		{
			return "ThresholdMower";
		}
  private:
    static const String info_;
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_THRESHOLDMOWER_H
