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
#ifndef OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
#define OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <vector>

namespace OpenMS
{
  /**
  	@brief Normalizer normalizes the peak intensities
 
	  @param method
	    0 = peaks get scaled relative to the maximum intensity<br>
	    1 = peaks get scaled relative to the TIC <br>
	    2 = peaks get scaled so that a vector with one dimension per peak has length 1<br>
	  @param windows
	    only applicable for method = 0 <br>
	    normalize in <i>windows</i> regions individually<br>
		
		@todo implement the method fully; implementing it correct; docs! (Andreas)
		@todo chenage SpectraFilter docu as soon as it is implemented fully (Andreas)
		
		@ingroup SpectraPreprocessing
  */
  class Normalizer
    : public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// standard constructor
    Normalizer();

    /// copy constructror
    Normalizer(const Normalizer& source);

    /// desctructor
    virtual ~Normalizer();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    Normalizer& operator = (const Normalizer& source);
		// @}

		// @name Accessors
		// @{
		///
    static PreprocessingFunctor* create() { return new Normalizer();}

		///
		static const String getName()
		{
			return "Normalizer";
		}

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
		
	    //std::vector<double> max = std::vector<double>((unsigned int)param_.getValue("windows"));
 	  	//double minmz = 10000;
    	//double maxmz = 0;
			
			////////////////////
			// just to normalize s.th.!! (normalize to one)
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

			//////////////////
			
			/*
    	for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
    	{
      	if ( it->getPosition()[0] < minmz ) minmz = it->getPosition()[0];
      	else if (it->getPosition()[0] > maxmz) maxmz = it->getPosition()[0];
    	}
    	for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
    	{
      	uint pos = (uint)((it->getPosition()[0] - minmz) / (maxmz - minmz) * max.size());
      	if (pos == max.size()) --pos;
      	if (max.at(pos) < it->getIntensity()) max.at(pos) = it->getIntensity();
    	}
    	for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
    	{
     		uint pos = (uint)((it->getPosition()[0] - minmz) / (maxmz - minmz) * max.size());
      	if (pos == max.size()) --pos;
      	it->setIntensity(it->getIntensity()/max[pos]);
    	}*/
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

    void filterPeakMap(PeakMap& exp);
		// @}

  };


}
#endif //OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
