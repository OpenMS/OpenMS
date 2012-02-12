// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_PARENTPEAKMOWER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <vector>

namespace OpenMS
{

  /**
  	@brief ParentPeakMower gets rid of high peaks that could stem from unfragmented precursor ions
		 
		@htmlinclude OpenMS_ParentPeakMower.parameters

		@ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI ParentPeakMower
		: public DefaultParamHandler 
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

		///
		template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
		{
			typedef typename SpectrumType::ConstIterator ConstIterator;
			typedef typename SpectrumType::Iterator Iterator;
			
			clean_all_charge_states_ = (Int)param_.getValue("clean_all_charge_states");
			consider_NH3_loss_ = (Int)param_.getValue("consider_NH3_loss");
			consider_H2O_loss_ = (Int)param_.getValue("consider_H2O_loss");
			window_size_ = (DoubleReal)param_.getValue("window_size");
			reduce_by_factor_ = (Int)param_.getValue("reduce_by_factor");
			factor_ = (DoubleReal)param_.getValue("factor");
			set_to_zero_ = (Int)param_.getValue("set_to_zero");

			if (spectrum.getMSLevel() == 1)
			{
				std::cerr << "Error: ParentPeakMower cannot be applied to MS level 1" << std::endl;
				return;
			}

    	//get precursor peak position precursorpeak
	    double pre_pos = 0.0;
			if (!spectrum.getPrecursors().empty()) pre_pos = spectrum.getPrecursors()[0].getMZ();
		
			if (pre_pos == 0)
			{
				std::cerr << "ParentPeakMower: Warning, Precursor Position not set" << std::endl;
				return;
			}

			Size pre_charge = spectrum.getPrecursors()[0].getCharge();
			if (pre_charge == 0)
			{
				default_charge_ = (Size)param_.getValue("default_charge");
				std::cerr << "ParentPeakMower: Warning, Precursor charge not set, assuming default charge (" << default_charge_ << ")" << std::endl;
				pre_charge = default_charge_;
			}

			pre_pos *= pre_charge;
					
			// identify the ranges which are to be considered
			std::vector<DRange<1> > ranges;
			for (Size z = 1; z <= pre_charge; ++z)
			{
				if (clean_all_charge_states_ || z == pre_charge)
				{
					// no adjusting needed for this charge
					DPosition<1> pre_z_pos, pos;
					DRange<1> range;
					
					// adjust the m/z by weight of precursor and charge
					pre_z_pos = DPosition<1>(pre_pos / double(z));
					range = DRange<1>(pre_z_pos - window_size_, pre_z_pos + window_size_);
					ranges.push_back(range);
										
					if (consider_NH3_loss_)
					{
						pos = DPosition<1>(pre_z_pos - 17.0 / double(z));
						range = DRange<1>(pos - window_size_, pos + window_size_);
						ranges.push_back(range);
					}
					if (consider_H2O_loss_)
					{
						pos = DPosition<1>(pre_z_pos - 18.0 / double(z));
						range = DRange<1>(pos - window_size_, pos + window_size_);
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
						if (reduce_by_factor_)
						{
							it->setIntensity(it->getIntensity() / factor_);
							break;
						}

						if (set_to_zero_)
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
		
		//TODO reimplement DefaultParamHandler::updateMembers_()

		//@}
		
	private:
		Size default_charge_;
		bool clean_all_charge_states_;
		bool consider_NH3_loss_;
		bool consider_H2O_loss_;
		DoubleReal window_size_;
		bool reduce_by_factor_;
		DoubleReal factor_;
		bool set_to_zero_;

  };

}
#endif // OPENMS_FILTERING/TRANSFORMERS_PARENTPEAKMOWER_H
