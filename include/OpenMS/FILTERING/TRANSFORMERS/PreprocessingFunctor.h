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

#ifndef OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H
#define OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
	/**
  	@brief Base class for Spectrum preprocessing classes
  	
  	@note Spectrum meta data arrays are invalidated by all preprocessing functors,
  	      that remove part of the peaks!
  */
  class OPENMS_DLLAPI PreprocessingFunctor
  	: public DefaultParamHandler
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
		
		/**
			@brief This method implements the actual filtering
			
			It is called by @ref filterSpectrum(SpectrumType&) and @ref filterPeakMap(PeakMap&) .
			It must be implemented in the derived classes!
			
			@note As it is a template method, this method is not usable through the base class pointer.
			      Through the base class pointer only @ref filterSpectrum(SpectrumType&) and @ref filterPeakMap(PeakMap&)
			      are usable. 
		*/
		template <typename SpectrumType> void filterSpectrum(SpectrumType& /*spectrum*/);

		/// filters an MSSpectrum, this method must be overwritten in the derived classes!
		virtual void filterPeakSpectrum(PeakSpectrum& spectrum) = 0;

		/// filters an MSExperiment, this method must be implemented in the derived classes!
		virtual void filterPeakMap(PeakMap& exp) = 0;
	};

}
#endif // OPENMS_FILTERING_TRANSFORMERS_PREPROCESSINGFUNCTOR_H
