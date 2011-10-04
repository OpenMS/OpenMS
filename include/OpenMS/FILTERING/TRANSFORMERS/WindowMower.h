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
#ifndef OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <set>

namespace OpenMS
{

  /**
  	@brief WindowMower augments the highest peaks in a sliding window
		 
		@htmlinclude OpenMS_WindowMower.parameters

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI WindowMower
		: public DefaultParamHandler 
  {
  public:

		// @name Constructors, destructors and assignment operators
		// @{
		/// default constructor
		WindowMower();
    /// destructor
    virtual ~WindowMower();

		/// copy constructor
		WindowMower(const WindowMower& source);
		/// assignment operator
		WindowMower& operator = (const WindowMower& source);
		// @}
		
		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
			
			windowsize_ = (DoubleReal)param_.getValue("windowsize");
    	peakcount_ = (UInt)param_.getValue("peakcount");

			//copy spectrum
			SpectrumType old_spectrum = spectrum;
			old_spectrum.sortByPosition();
			
			//find high peak positions
			bool end  = false;
			std::set<double> positions;
			for (ConstIterator it = old_spectrum.begin(); it != old_spectrum.end(); ++it)
			{
				// copy the window from the spectrum
				SpectrumType window;
				for (ConstIterator it2 = it; (it2->getPosition() - it->getPosition() < windowsize_); )
				{
					window.push_back(*it2);
					if (++it2 == old_spectrum.end())
					{
						end = true;
						break;
					}
				}
				
				//extract peakcount most intense peaks				
				window.sortByIntensity(true);
				for (Size i = 0; i < peakcount_; ++i)
				{
					if (i < window.size())
					{
						positions.insert(window[i].getMZ());
					}
				}
				//abort at the end of the spectrum
				if (end) break;
			}

			// replace the old peaks by the new ones
			spectrum.clear(false);
			for (ConstIterator it = old_spectrum.begin(); it != old_spectrum.end(); ++it)
			{
				if (positions.find(it->getMZ()) != positions.end())
				{
					spectrum.push_back(*it);
				}
			}
    }

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);

		//TODO reimplement DefaultParamHandler::updateMembers_()
		
		private:
			DoubleReal windowsize_;
			UInt peakcount_;
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
