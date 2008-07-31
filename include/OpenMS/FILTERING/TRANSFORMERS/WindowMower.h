// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <set>

namespace OpenMS
{

  /**
  	@brief WindowMower augments the highest peaks in a sliding window
		 
		@ref WindowMower_Parameters are explained on a separate page.

		@ingroup SpectraPreprocessers
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
    static PreprocessingFunctor* create() { return new WindowMower(); }
		
		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::Iterator Iterator;
			typedef typename SpectrumType::ConstIterator ConstIterator;
			typedef typename SpectrumType::ContainerType ContainerType;
			
			double windowsize = (double)param_.getValue("windowsize");
    	UInt peakcount = (int)param_.getValue("peakcount");

			std::set<double> positions; // store the indices that are the most intense ones in an interval
			
			spectrum.getContainer().sortByPosition();
			
			// slide the window over spectrum and store peakcount most intense peaks (if available) of every window position
			bool end(false);
			for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
			{
				ContainerType container;
				for (ConstIterator it2 = it; (it2->getPosition() - it->getPosition() < windowsize); )
				{
					container.push_back(*it2);
					if (++it2 == spectrum.end())
					{
						end = true;
						break;
					}
				}

				container.sortByIntensity(true);
				
				for (UInt i = 0; i < peakcount; ++i)
				{
					if (container.size() > i)
					{
						positions.insert(container[i].getMZ());
					}
				}

				if (end)
				{
					break;
				}
			}

			// add the found peaks to a new container
			ContainerType container;
			for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
			{
				if (positions.find(it->getMZ()) != positions.end())
				{
					container.push_back(*it);
				}
			}
			
			// overwrite the spectrum with the new container
			spectrum.setContainer(container);
			spectrum.getContainer().sortByPosition();
    }

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);
	
		///
		static const String getProductName()
		{
			return "WindowMower";
		}
		// @}
		
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
