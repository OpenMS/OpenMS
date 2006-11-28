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
  	@brief ParentPeakMower gets rid of high peaks that could stem from unfragmented precursor ions
  
	  @param windowsize consider all peaks inside parent ion m/z +- windowsize

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
    static PreprocessingFunctor* create() { return new ParentPeakMower(); }

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
    	//get precursor peak position precursorpeak
    	double pppos = spectrum.getPrecursorPeak().getPosition()[0];
			//???? Why devide by charge? Scaling is done in the wrong place like that?!?!?!
			if (spectrum.getPrecursorPeak().getCharge() != 0)
			{
				pppos /= spectrum.getPrecursorPeak().getCharge();
			}
			//abort if no precursor peak was set
			if (pppos==0.0)
			{
				if (spectrum.getMSLevel()>1)
				{
					std::cout << "Warning: No precursor peak for Spectrum with MS-level greater than '1'. Aborting!" << std::endl;
				}
				return;
			}

    	spectrum.getContainer().sortByPosition();

    	// calculate mean
    	double mean = 0;
    	for (typename SpectrumType::ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
    	{
      	mean += it->getIntensity();
    	}
    	mean /= spectrum.size();
			
			//do scaling
			double window = (double)param_.getValue("windowsize");
    	for (typename SpectrumType::Iterator it = spectrum.MZBegin(pppos - window); it != spectrum.MZEnd(pppos + window); ++it)
    	{
				if (it->getIntensity() > mean) it->setIntensity(mean);
    	}
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);

		///
		static const String getName()
		{
			return "ParentPeakMower";
		}
		//@}
  };

}
#endif // OPENMS_FILTERING/TRANSFORMERS_PARENTPEAKMOWER_H
