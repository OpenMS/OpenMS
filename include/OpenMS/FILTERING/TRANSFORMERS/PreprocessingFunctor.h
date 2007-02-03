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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H
#define OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
	/**
  	@brief Base class for Spectrum preprocessing classes
  */
  class PreprocessingFunctor : public FactoryProduct
  {
  public:

    /// default constructor
    PreprocessingFunctor();
		
    /// copy constructor
    PreprocessingFunctor(const PreprocessingFunctor& source);
		
    /// destructor
    virtual ~PreprocessingFunctor();
		
    /// assignment operator
    PreprocessingFunctor& operator = (const PreprocessingFunctor& source);

		static void registerChildren();
		
		/// this is just an interface method, it must be implemented in the derived classes
		template <typename SpectrumType> void filterSpectrum(SpectrumType& /*spectrum*/);

		/// filters an MSSpectrum, this method should be overwritten in the derived classes
		virtual void filterPeakSpectrum(PeakSpectrum& spectrum) = 0;

		/// filters an MSExperiment, this method should be overwritten in the derived classes
		virtual void filterPeakMap(PeakMap& exp) = 0;

		///
		static const String getProductName()
		{
			return "PreprocessingFunctor";
		}

	};

}
#endif // OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H
