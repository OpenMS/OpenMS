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
#ifndef OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#include <map>
#include <vector>
#include <cmath>
#include <utility>

namespace OpenMS
{

  /**
	  @brief IsotopeMarker marks peak pairs which could represent an ion and its isotope

		@todo implement a real isotope marking here with isotopedistributions and fitting (Andreas)

		@htmlinclude OpenMS_IsotopeMarker.parameters

		@ingroup PeakMarker
  */
  class OPENMS_DLLAPI IsotopeMarker
    : public PeakMarker
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    IsotopeMarker();

    /// copy constructor
    IsotopeMarker(const IsotopeMarker& source);

    /// destructor
    virtual ~IsotopeMarker();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    IsotopeMarker& operator=(const IsotopeMarker& source);
		// @}

		// @name Accessors
		// @{
		///
    static PeakMarker* create() { return new IsotopeMarker();}

		///
		template <typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum)
		{
			double mzvariation = (double)param_.getValue("mz_variation");
    	double invariation = (double)param_.getValue("in_variation");
    	Size marks = param_.getValue("marks");
			
    	spectrum.sortByPosition();
			
			std::map<double, Size> isotopemarks ; // possible isotopes
			
    	for (Size i = 0; i < spectrum.size(); ++i)
    	{
      	double mz = spectrum[i].getPosition()[0];
      	double intensity = spectrum[i].getIntensity();
      	Size j = i+1;

      	//std::vector<std::pair<double, double> > isotopes = SpectrumGenerator::instance()->isotopepeaks(mz, intensity);
				IsotopeDistribution id;
				id.estimateFromPeptideWeight(mz);

      	while (j < spectrum.size() && spectrum[j].getPosition()[0] <= mz + 3 + mzvariation)
      	{
        	double curmz = spectrum[j].getPosition()[0];
        	double curIntensity = spectrum[j].getIntensity();
        	UInt iso = (UInt)(curmz - mz + 0.499999);
        	if (iso > 0 && curmz - mz - iso > mzvariation)
        	{
          	++j;
          	continue;
        	}
        	if (std::fabs(id.begin()->second * intensity - curIntensity) < invariation * id.begin()->second * intensity)
        	{
          	isotopemarks[mz]++;
        	  isotopemarks[curmz]++;
        	}
        	++j;
      	}
    	}
			
    	for (std::map<double, Size>::const_iterator cmit = isotopemarks.begin(); cmit != isotopemarks.end(); ++cmit)
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
			return "IsotopeMarker";
		}
		// @}
		
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H
