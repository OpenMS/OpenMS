// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <vector>

namespace OpenMS
{

  /**
  	@brief ParentPeakMower gets rid of high peaks that could stem from unfragmented precursor ions
		 
		@htmlinclude OpenMS_ParentPeakMower.parameters

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI ParentPeakMower : public PreprocessingFunctor
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
			typedef typename SpectrumType::ConstIterator ConstIterator;
			typedef typename SpectrumType::Iterator Iterator;
			
			if (spectrum.getMSLevel() == 1)
			{
				std::cerr << "Error: ParenPeakMower cannot be applied to MS-level 1" << std::endl;
				return;
			}

    	//get precursor peak position precursorpeak
	    double pre_pos = 0.0;
			if (!spectrum.getPrecursors().empty()) pre_pos = spectrum.getPrecursors()[0].getMZ();
		
			if (pre_pos == 0)
			{
				std::cerr << "ParenPeakMower: Warning, Precursor Position not set" << std::endl;
				return;
			}

			UInt pre_charge = spectrum.getPrecursors()[0].getCharge();
			if (pre_charge == 0)
			{
				UInt default_charge = (unsigned int)param_.getValue("default_charge");
				std::cerr << "ParentPeakMower: Warning, Precursor charge not set, assuming default charge (" << default_charge << ")" << std::endl;
				pre_charge = default_charge;
			}

			pre_pos *= pre_charge;
			
			// get all other parameters 
			bool clean_all_charge_states = (Int)param_.getValue("clean_all_charge_states");
			bool consider_NH3_loss = (Int)param_.getValue("consider_NH3_loss");
			bool consider_H2O_loss = (Int)param_.getValue("consider_H2O_loss");
			double window_size = (double)param_.getValue("window_size");
			bool reduce_by_factor = (Int)param_.getValue("reduce_by_factor");
			double factor = (double)param_.getValue("factor");
			bool set_to_zero = (Int)param_.getValue("set_to_zero");
		
			// identify the ranges which are to be considered
			std::vector<DRange<1> > ranges;
			for (Size z = 1; z <= pre_charge; ++z)
			{
				if (clean_all_charge_states || z == pre_charge)
				{
					// no adjusting needed for this charge
					DPosition<1> pre_z_pos, pos;
					DRange<1> range;
					
					// adjust the m/z by weight of precursor and charge
					pre_z_pos = DPosition<1>(pre_pos / double(z));
					range = DRange<1>(pre_z_pos - window_size, pre_z_pos + window_size);
					ranges.push_back(range);
										
					if (consider_NH3_loss)
					{
						pos = DPosition<1>(pre_z_pos - 17.0 / double(z));
						range = DRange<1>(pos - window_size, pos + window_size);
						ranges.push_back(range);
					}
					if (consider_H2O_loss)
					{
						pos = DPosition<1>(pre_z_pos - 18.0 / double(z));
						range = DRange<1>(pos - window_size, pos + window_size);
						ranges.push_back(range);
					}
				}
			}

//			for (std::vector<DRange<1> >::const_iterator rit = ranges.begin(); rit != ranges.end(); ++rit)
//			{
//				std::cerr << *rit << std::endl;
//			}
			
			// apply the intensity reduction to the collected ranges
			for (Iterator it = spectrum.begin(); it != spectrum.end(); ++it)
			{
				for (std::vector<DRange<1> >::const_iterator rit = ranges.begin(); rit != ranges.end(); ++rit)
				{
					if (rit->encloses(it->getPosition()))
					{
						if (reduce_by_factor)
						{
							it->setIntensity(it->getIntensity() / factor);
							break;
						}

						if (set_to_zero)
						{
							it->setIntensity(0.0);
							break;
						}
					}
				}
			}
			
			return;
		}

		void filterPeakSpectrum(PeakSpectrum& spectrum);

		void filterPeakMap(PeakMap& exp);

		///
		static const String getProductName()
		{
			return "ParentPeakMower";
		}
		//@}
  };

}
#endif // OPENMS_FILTERING/TRANSFORMERS_PARENTPEAKMOWER_H
