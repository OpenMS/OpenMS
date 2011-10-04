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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <cmath>

namespace OpenMS
{
  /**
  	@brief Scales the intensity of peaks to the sqrt 

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI SqrtMower
 		: public DefaultParamHandler 
 {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    SqrtMower();
    /// destructor
    virtual ~SqrtMower();
	
		/// copy constructor
		SqrtMower(const SqrtMower& source);
		/// assignment operator
		SqrtMower& operator=(const SqrtMower& source);
		// @}

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{      
      bool warning = false;
			for (typename SpectrumType::Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        DoubleReal intens = it->getIntensity();        
        if(intens<0)
        {
          intens=0;
          warning=true;
        }
        it->setIntensity(std::sqrt(intens));
			}
      if(warning)
      {
        std::cerr<<"Warning negative intensities were set to zero"<<std::endl;
      }
			return;
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);
		
		//TODO reimplement DefaultParamHandler::updateMembers_()
		
		// @}
		
  };

}

#endif // OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H
