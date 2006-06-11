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
// $Id: Normalizer.h,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
#define OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H

#include <OpenMS/FILTERING/TRANSFORMERS/MowerFunctor.h>
#include <vector>

namespace OpenMS
{
  /**
  	@brief Normalizer normalizes the peak intensities
 
		@todo implement the method fully; implementing it correct; docs!
 
	  \param method
	    0 = peaks get scaled relative to the maximum intensity<br>
	    1 = peaks get scaled relative to the TIC <br>
	    2 = peaks get scaled so that a vector with one dimension per peak has length 1<br>
	  \param windows
	    only applicable for method = 0 <br>
	    normalize in <i>windows</i> regions individually<br>
  */
  class Normalizer
    :public MowerFunctor
  {
  public:
    /// standard constructor
    Normalizer();

    /// copy constructror
    Normalizer(const Normalizer& source);

    /// desctructor
    ~Normalizer();

    /// assignment operator
    Normalizer& operator=(const Normalizer& source);

    static FactoryProduct* create() { return new Normalizer();}
    //void operator()(MSSpectrum< DPeak<1> >&) const;
    //String info() const;

		static const String getName()
		{
			return "Normalizer";
		}


		template <typename SpectrumType> void apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
		
	    std::vector<double> max = std::vector<double>((unsigned int)param_.getValue("windows"));
 	  	double minmz = 10000;
    	double maxmz = 0;
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
    	}
		}

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_NORMALIZER_H
