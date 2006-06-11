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
// $Id: MowerFunctor.h,v 1.7 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_MOWERFUNCTOR_H
#define OPENMS_FILTERING_TRANSFORMERS_MOWERFUNCTOR_H

#include <vector>
#include <map>
#include <string>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{
  /**
  MowerFunctors do some form of preprocessing to a MSSpectrum<br>
  todo: rename to PreprocessFunctor
  */
  class MowerFunctor : public FactoryProduct
  {
  public:
    /** @brief standard constructor <br> */
    MowerFunctor() ;
    /** @brief copy constructor <br> */
    MowerFunctor(const MowerFunctor& source);
    /** @brief destructor <br> */
    virtual ~MowerFunctor() {}
    /** @brief assignment operator < br> */
    MowerFunctor& operator=(const MowerFunctor& source);

    //virtual void operator()(Spectrum< DPeak<1> >&) const = 0;

		template <typename SpectrumType> void operator()(SpectrumType&) {}
	};

}
#endif // OPENMS_FILTERING_TRANSFORMERS_MOWERFUNCTOR_H
