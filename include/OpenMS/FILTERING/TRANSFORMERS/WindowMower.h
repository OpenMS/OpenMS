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
// $Id: WindowMower.h,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <map>

namespace OpenMS
{

  /**
  	@brief WindowMower augments the highest peaks in a sliding window <br>
  
	  \param peakcount: nr of peaks that are augmented in each step
	  \param windowsize: size of sliding window
  */
  class WindowMower
  	 : public PreprocessingFunctor
  {
  public:
    /// standard constructor
    WindowMower();

    /// copy constructor
    WindowMower(const WindowMower& source);

    /// destructor
    ~WindowMower();

    /// assignment operator
    WindowMower& operator=(const WindowMower& source);

    static FactoryProduct* create() { return new WindowMower();}
		
    //void operator()(Spectrum< DPeak<1> >&) const;
	
		template <typename SpectrumType> void apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
			
			double windowsize = (double)param_.getValue("windowsize");
    	uint peakcount = (int)param_.getValue("peakcount");
    	std::map<double, double> peaksinwindow; // peakheight,pos
    	std::map<double, int> marks; // peaks get marked if they belong to the <peakcount> highest in the window
			
			for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
    	{
      	peaksinwindow.clear();
       	for (uint i = 0; (it+i) != spectrum.end() && (it+i)->getIntensity() < it->getIntensity()+windowsize ; ++i)
       	{
        	 peaksinwindow.insert(std::make_pair<double, double>((it+i)->getIntensity(), (it+i)->getPosition()[0]));
       	}
			 
       	std::map<double, double>::reverse_iterator it2 = peaksinwindow.rbegin();
       	for (uint i = 0; i < peakcount && i < peaksinwindow.size(); ++i)
       	{
         	marks[(it2++)->second]++;
       	}
			 
       	//todo do something with the marking, maybe multiply the peaks with it
       	for (uint i = 0; i < spectrum.size(); ++i)
       	{
         	spectrum.getContainer()[i].setIntensity(spectrum.getContainer()[i].getIntensity() + spectrum.getContainer()[i].getIntensity()*marks[spectrum.getContainer()[i].getPosition()[0]]);
       	}
      }    
    }
	
    //String info() const;
    void operator()(MSSpectrum< DPeak<1> >&) const;
    String info() const;

		static const String getName()
		{
			return "WindowMower";
		}
  private:
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
