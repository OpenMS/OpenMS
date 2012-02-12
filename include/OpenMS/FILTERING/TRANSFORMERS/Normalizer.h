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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
#define OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <vector>

namespace OpenMS
{
  /**
  	@brief Normalizer normalizes the peak intensities
		 
		@htmlinclude OpenMS_Normalizer.parameters
		
		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI Normalizer
		: public DefaultParamHandler 
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    Normalizer();
    /// desctructor
    virtual ~Normalizer();
	
		/// assignment operator
    Normalizer& operator = (const Normalizer& source);
    /// copy constructror
    Normalizer(const Normalizer& source);
		
		// @}

		// @name Accessors
		// @{

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
				
			method_ = param_.getValue("method");

			// normalizes the max peak to 1 and the rest of the peaks to values relative to max
			if (method_ == "to_one")
			{
				double max(0);
				for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
				{
					if (max < it->getIntensity())
					{
						max = it->getIntensity();
					}
				}
				for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
				{
					it->setIntensity(it->getIntensity() / max);
				}
			}
			// normalizes the peak intensities to the TIC
			else if (method_ == "to_TIC")
			{
				double sum(0);
				for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
				{
					sum += it->getIntensity();
				}

				for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
				{
					it->setIntensity(it->getIntensity() / sum);
				}
			}
 			// method unknown
      else
      {
  			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Method not known", method_);
      }
			return;
			
		}

		///
		void filterPeakSpectrum(PeakSpectrum& spectrum);
		///
    void filterPeakMap(PeakMap& exp);

		//TODO reimplement DefaultParamHandler::updateMembers_()

		// @}
		
	private:
		String method_;
	
  };


}
#endif //OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
