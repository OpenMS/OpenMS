// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

namespace OpenMS
{

  /**
  	@brief ParentPeakMower gets rid of high peaks that could stem from unfragmented parent ions<br>
  
	  \param windowsize consider all peaks inside parent ion m/z +- windowsize
	  \param x what is considered high: intensity > x*mean(peakintensity)

		@ingroup SpectraPreprocessing
  */
  class ParentPeakMower : public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    ParentPeakMower();

    /// copy constructor
    ParentPeakMower(const ParentPeakMower& source);

    /// destructor
    virtual ~ParentPeakMower();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    ParentPeakMower& operator = (const ParentPeakMower& source);
		// @}

		// @name Accessors
		// @{
		///
    static FactoryProduct* create() { return new ParentPeakMower(); }

		///
		template <typename SpectrumType> void apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::ConstIterator ConstIterator;
			typedef typename SpectrumType::Iterator Iterator;
		
    	double window = (double)param_.getValue("windowsize");
    	double mean = 0;

    	spectrum.getContainer().sortByPosition();

    	// calculate mean
    	for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
    	{
      	mean += it->getIntensity();
    	}
    	mean /= spectrum.size();

    	// assumed position of precursorpeak
    	double pppos = spectrum.getPrecursorPeak().getPosition()[0];
			if (spectrum.getPrecursorPeak().getCharge() != 0)
			{
				pppos /= spectrum.getPrecursorPeak().getCharge();
			}

    	Iterator it = spectrum.end();
			
    	if (it == spectrum.begin()) return;
    	do
    	{
      	--it;
      	if (it->getPosition()[0] <= pppos + window && it->getPosition()[0] >= pppos - window)
      	{
        	if (it->getIntensity() > mean) it->setIntensity(mean);
      	}
      	else if (it->getPosition()[0] < pppos - window)
      	{
        	break;
      	}
    	}	while (it != spectrum.begin());
		}

		///
		static const String getName()
		{
			return "ParentPeakMower";
		}
		//@}
  };

}
#endif // OPENMS_FILTERING/TRANSFORMERS_PARENTPEAKMOWER_H
