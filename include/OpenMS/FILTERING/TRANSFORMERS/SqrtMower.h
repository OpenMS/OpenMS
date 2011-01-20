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
#ifndef OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <cmath>

namespace OpenMS
{
  /**
  	@brief Scales the intensity of peaks to the sqrt 

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI SqrtMower : public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    SqrtMower();

    /// copy constructor
    SqrtMower(const SqrtMower& source);

    /// destructor
    virtual ~SqrtMower();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    SqrtMower& operator=(const SqrtMower& source);
		// @}

		// @name Accessors
		// @{
		///
    static PreprocessingFunctor* create() { return new SqrtMower(); }

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			for (typename SpectrumType::Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
			{
				it->setIntensity(std::sqrt(it->getIntensity()));
			}
			return;
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);

		///
		static const String getProductName()
		{
			return "SqrtMower";
		}
		// @}
		
  };

}

#endif // OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H
