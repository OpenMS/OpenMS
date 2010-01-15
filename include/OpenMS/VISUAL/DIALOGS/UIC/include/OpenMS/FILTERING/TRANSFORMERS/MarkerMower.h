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
#ifndef OPENMS_FILTERING_TRANSFORMERS_MARKERMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_MARKERMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <vector>
#include <map>

namespace OpenMS
{
  /**
  	@brief MarkerMower uses PeakMarker to find peaks, those that are not marked get removed

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI MarkerMower : public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    MarkerMower();

    /// copy constructor
    MarkerMower(const MarkerMower& source);

    /// destructor
    virtual ~MarkerMower();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    MarkerMower& operator = (const MarkerMower& source);
		// @}

		// @name Accessors
		// @{
		///
    static PreprocessingFunctor* create() { return new MarkerMower(); }

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
		
			std::map<double, int> marks;
    	for (std::vector<PeakMarker*>::const_iterator cvit = markers_.begin(); cvit != markers_.end(); ++cvit)
    	{
      	std::map<double, bool> marked;
				(*cvit)->apply(marked, spectrum);
      	for (std::map<double, bool>::const_iterator cmit = marked.begin(); cmit != marked.end(); ++cmit)
      	{
        	if (cmit->second) 
					{
						marks[cmit->first]++;
					}
      	}
    	}

			for (Iterator it = spectrum.begin(); it != spectrum.end(); )
			{
 				if (marks[it->getMZ()] > 0)
				{
					++it;
				}
				else
				{
					it = spectrum.erase(it);
				}
			}
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);
		
		static const String getProductName()
		{
			return "MarkerMower";
		}

    /// insert new Marker (violates the PreprocessingFunctor interface)
    void insertmarker(PeakMarker* peak_marker);
		// @}
	
	private: 
	
    /// used peak markers
    std::vector<PeakMarker*> markers_;
  };
}
#endif // OPENMS_COMPARISON_CLUSTERING_MARKERMOWER_H
