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
#ifndef OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <vector>
#include <algorithm>
#include <set>

namespace OpenMS
{

  /**
  	@brief WindowMower augments the highest peaks in a sliding window <br>
  
	  @param peakcount: nr of peaks that are augmented in each step
	  @param windowsize: size of sliding window

		@ingroup SpectraPreprocessing
  */
  class WindowMower
  	 : public PreprocessingFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    WindowMower();

    /// copy constructor
    WindowMower(const WindowMower& source);

    /// destructor
    virtual ~WindowMower();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    WindowMower& operator = (const WindowMower& source);
		// @}

		// @name Accessors
		// @{
		///
    static FactoryProduct* create() { return new WindowMower(); }
		
		///
		template <typename SpectrumType> void apply(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
			typedef typename SpectrumType::ContainerType ContainerType;
			
			double windowsize = (double)param_.getValue("windowsize");
    	Size peakcount = (int)param_.getValue("peakcount");

			std::set<double> positions; // store the indices that are the most intense ones in an interval
			
			spectrum.getContainer().sortByPosition();
			
			// slide the window over spectrum and store peakcount most intense peaks (if available) of every window position
			for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
			{
				ContainerType container;
				for (Iterator it2 = it; (it2->getPosition()[0] - it->getPosition()[0] < windowsize) && it2 != spectrum.end(); ++it2)
				{
					container.push_back(*it2);
				}

				container.sortByIntensity();

				for (Size i = 0; i != peakcount; ++i)
				{
					if (container.size() > i)
					{
						positions.insert(container[i].getPosition()[0]);
					}
				}
			}

			// add the found peaks to a new container
			ContainerType container;
			for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
			{
				if (positions.find(it->getPosition()[0]) != positions.end())
				{
					container.push_back(*it);
				}
			}
			
			// overwrite the spectrum with the new container
			spectrum.setContainer(container);

			/* old implementation ?!?
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
      }*/
    }
	
		///
		static const String getName()
		{
			return "WindowMower";
		}
		// @}
		
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
