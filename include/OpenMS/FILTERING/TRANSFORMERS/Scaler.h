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
#ifndef OPENMS_FILTERING_TRANSFORMERS_SCALER_H
#define OPENMS_FILTERING_TRANSFORMERS_SCALER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <map>

namespace OpenMS
{
  /**
  	@brief Scaler scales the peak by ranking the peaks and assigning intensity according to rank<br>

		@ingroup SpectraPreprocessing
  */
  class Scaler
    : public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    Scaler();

    /// copy constructor
    Scaler(const Scaler& source);

    /// destructor
    virtual ~Scaler();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    Scaler& operator = (const Scaler& source);
		// @}

		// @name Accessors
		// @{
		///
    static FactoryProduct* create() { return new Scaler(); }

		///
		static const String getName()
		{
			return "Scaler";
		}

		template <typename SpectrumType> void apply(SpectrumType& spectrum)
		{	
			typedef typename SpectrumType::ConstIterator ConstIterator;
			typedef typename SpectrumType::Iterator Iterator;
			std::map<double, int> peakssorted;
	    int count(0);
			
 	  	for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
 	   	{
 	    	peakssorted[it->getIntensity()] = 0;
 	    	++count;
 	  	}
			
    	for(std::map<double, int>::reverse_iterator rmit = peakssorted.rbegin(); rmit != peakssorted.rend(); ++rmit)
    	{
      	if (--count > 0) 
				{
					peakssorted[rmit->first] = count;
				}
    	}
			
    	for (Iterator it = spectrum.begin(); it != spectrum.end(); )
    	{
      	double new_intensity = peakssorted[it->getIntensity()];
      	if (new_intensity > 0)
      	{
        	it->setIntensity(peakssorted[it->getIntensity()]);
        	++it;
      	}
      	else
      	{
        	it = spectrum.getContainer().erase(it);
      	}
    	}
		}
		// @}

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_SCALER_H
