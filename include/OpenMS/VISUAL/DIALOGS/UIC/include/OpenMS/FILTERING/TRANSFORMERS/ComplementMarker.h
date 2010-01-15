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
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief ComplementMarker marks peak pairs which could represent y - b ion pairs
    
		@htmlinclude OpenMS_ComplementMarker.parameters

		@ingroup PeakMarker
  */
  class OPENMS_DLLAPI ComplementMarker
    : public PeakMarker
  {
  public:

		// @name Constructors and Destructors
		//@{
    /// standard constructor
    ComplementMarker();

    /// copy constructor
    ComplementMarker(const ComplementMarker& source);

    /// destructor
    virtual ~ComplementMarker();
		//@}

		// @name Operators
		//@{
    /// assignment operator
    ComplementMarker& operator = (const ComplementMarker& source);
		//@}

		// @name Accessors
		//@{
		///
    static PeakMarker* create() { return new ComplementMarker(); }
		
		///
		template <typename SpectrumType> void apply(std::map<double, bool> marked, SpectrumType& spectrum)
		{
			if (spectrum.size() < 2)
			{
				return;
			}
		
			// how often a peak needs to be marked to be returned
    	double marks = (double)param_.getValue("marks");
	    double parentmass = 0.0;
			if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
    	double tolerance = (double)param_.getValue("tolerance");
    	std::map<double, int> matching_b_y_ions;
			
    	spectrum.sortByPosition();
			
    	SignedSize j = spectrum.size() -1;
    	for (Size i = 0; i < spectrum.size(); ++i)
    	{
      	while (j >= 0 && spectrum[j].getPosition()[0] > (parentmass - spectrum[i].getPosition()[0]) + tolerance)
				{
        	j--;
      	}
				
      	// just takes the first matching ion; todo take all
      	if (j >= 0 && std::fabs(spectrum[i].getPosition()[0] + spectrum[j].getPosition()[0] - parentmass) < tolerance)
      	{
        	matching_b_y_ions[spectrum[i].getPosition()[0]]++;
        	matching_b_y_ions[spectrum[j].getPosition()[0]]++;
        	j--;
      	}
    	}

			for (std::map<double, int>::const_iterator cmit = matching_b_y_ions.begin(); cmit != matching_b_y_ions.end(); ++cmit)
			{
				if (cmit->second >= marks)
				{
					marked.insert(std::make_pair<double, bool>(cmit->first, true));
				}
			}
		}
	
		/// returns the name to register at the factory
		static const String getProductName()
		{
			return "ComplementMarker";
		}
		//@}
		
  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTMARKER_H
