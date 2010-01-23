// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
  	@brief NeutralLossMarker marks peak pairs which could represent an ion an its neutral loss (water, ammonia)
		 
		@htmlinclude OpenMS_NeutralLossMarker.parameters

		@ingroup PeakMarker
  */
  class OPENMS_DLLAPI NeutralLossMarker
    : public PeakMarker
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    NeutralLossMarker();

    /// copy constructor
    NeutralLossMarker(const NeutralLossMarker& source);

    /// destructor
    virtual ~NeutralLossMarker();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    NeutralLossMarker& operator = (const NeutralLossMarker& source);
		// @}

		// @name Accessors
		// @{
		///
    static PeakMarker* create() { return new NeutralLossMarker(); }

		///
		template <typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum)
		{
			// how often a peak needs to be marked to be returned
    	double marks = (double)param_.getValue("marks");
    	double tolerance = (double)param_.getValue("tolerance");
    	std::map<double, SignedSize> ions_w_neutrallosses;
    	spectrum.sortByPosition();
    	for (Size i = 0; i < spectrum.size(); ++i)
    	{
      	double mz = spectrum[i].getPosition()[0];
      	double intensity = spectrum[i].getIntensity();
      	SignedSize j = i - 1;
      	while (j >= 0)
      	{
        	double curmz = spectrum[j].getPosition()[0];
        	double curIntensity = spectrum[j].getIntensity();

        	// check for peak thats a a water or ammonia away
        	if (std::fabs(mz - curmz - 17) < tolerance || std::fabs(mz - curmz - 18) < tolerance)
        	{
          	// neutral loss peak should be smaller
          	if (curIntensity < intensity)
          	{
            	ions_w_neutrallosses[mz]++;
            	// neutral loss peak not marked
            	//ions_w_neutrallosses[curmz]++;
          	}
        	}
        	else 
					{
						if (mz - curmz > 18.3)
		        {
    		      break;
        		}
					}
        	--j;
      	}
    	}
			
    	for (std::map<double, SignedSize>::const_iterator cmit = ions_w_neutrallosses.begin(); cmit != ions_w_neutrallosses.end(); ++cmit)
    	{
      	if (cmit->second >= marks)
      	{
        	marked.insert(std::make_pair<double, bool>(cmit->first, true));
      	}
    	}
    	return;
		}

		///
		static const String getProductName()
		{
			return "NeutralLossMarker";
		}
		// @}

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H
