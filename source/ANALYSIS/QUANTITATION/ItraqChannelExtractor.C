// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/KERNEL/RangeUtils.h>

namespace OpenMS
{
	/// Constructor (assuming 4plex)
	ItraqChannelExtractor::ItraqChannelExtractor()
		:DefaultParamHandler("ItraqChannelExtractor"),
		 itraq_type_(ItraqConstants::FOURPLEX),
		 channel_map_()
	{
		init_();
		setDefaultParams_();
	}

	/// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES)
	ItraqChannelExtractor::ItraqChannelExtractor(Int itraq_type)
		:DefaultParamHandler("ItraqChannelExtractor"),
		 itraq_type_(itraq_type),
		 channel_map_()
	{
		init_();
		setDefaultParams_();
	}

	/// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES) and param
	ItraqChannelExtractor::ItraqChannelExtractor(Int itraq_type, const Param& param)
		:DefaultParamHandler("ItraqChannelExtractor"),
		 itraq_type_(itraq_type),
		 channel_map_()
	{
		init_();
		setDefaultParams_();
		setParameters(param);
		updateMembers_();
	}

	/// copy constructor
	ItraqChannelExtractor::ItraqChannelExtractor(const ItraqChannelExtractor& cp)
	: DefaultParamHandler(cp),
		ItraqConstants(cp),
		itraq_type_(cp.itraq_type_), 
		channel_map_(cp.channel_map_)
	{

	}

	/// assignment operator
	ItraqChannelExtractor& ItraqChannelExtractor::operator = (const ItraqChannelExtractor& rhs)
	{
		if (this == &rhs) return *this;

		DefaultParamHandler::operator = (rhs);
		ItraqConstants::operator = (rhs);
		itraq_type_ = rhs.itraq_type_; 
		channel_map_ = rhs.channel_map_;

		return *this;
	}


	/// @brief extracts the iTRAQ channels from the MS data and stores intensity values in a consensus map
	///
	/// @param ms_exp_data Raw data to read
	/// @param consensus_map Output each MS² scan as a consensus feature
	void ItraqChannelExtractor::run(const MSExperiment<Peak1D>& ms_exp_data, ConsensusMap& consensus_map)
	{
		MSExperiment<> ms_exp_MS2;

		//check for peak type (raw data required)
		bool do_picking = false;
		if ((String) param_.getValue("intensity_method")=="peak_picker")
		{
			do_picking = true;
			if (ms_exp_data.size()>0)
			{
				for (Size i=0; i<ms_exp_data[0].getDataProcessing().size(); ++i)
				{
					if (ms_exp_data[0].getDataProcessing()[i].getProcessingActions().count(DataProcessing::PEAK_PICKING)==1)
					{
						std::cerr << "WARNING: ItraqChannelExtractor::run() MSExperiment is already picked! Skipping PeakPicking!" << std::endl;
						do_picking = false;
						break;
					}
				}
			}
		}
		else
		{
			do_picking = false;
		}

		if (do_picking)
		{
			#ifdef ITRAQ_DEBUG
			std::cout << "picking peaks..." << std::endl;
			#endif
			// pick
			if (PeakTypeEstimator().estimateType(ms_exp_data[0].begin(),ms_exp_data[0].end())==SpectrumSettings::PEAKS)
			{
				std::cerr << "Warning: OpenMS peak type estimation indicates that this is peak data already!";
			}					

			Param pepi_param = param_.copy("PeakPicker:",true);
			PeakPickerCWT peak_picker;
			peak_picker.setParameters(pepi_param);

			for (size_t idx=0; idx<ms_exp_data.size(); ++idx)
			{
				if (ms_exp_data[idx].getMSLevel()==2)
				{
					// copying is important, because the peakpicker does not preserve precusor peak information etc
					MSSpectrum <Peak1D> spectrum_picked(ms_exp_data[idx]);
					peak_picker.pick(ms_exp_data[idx],spectrum_picked);
					ms_exp_MS2.push_back(spectrum_picked);
				}
			}
		}
		else
		{ // copy only MS² scans
			for (size_t idx=0; idx<ms_exp_data.size(); ++idx)
			{
				if (ms_exp_data[idx].getMSLevel()==2)
				{
					ms_exp_MS2.push_back(ms_exp_data[idx]);
				}
			}				
		}
		#ifdef ITRAQ_DEBUG
		std::cout << "we have " << ms_exp_MS2.size() << " scans left of level " << ms_exp_MS2[0].getMSLevel() << std::endl;
		std::cout << "run: channel_map_ has " << channel_map_.size() << " entries!" << std::endl;
		#endif
		consensus_map.clear(false);
		// set <mapList> header
		Int index = 0;
		for (ChannelMapType::const_iterator cm_it = channel_map_.begin(); cm_it!=channel_map_.end(); ++cm_it)
		{
			// structure of Map cm_it
			//   first == channel-name as Int e.g. 114
			//	 second == ChannelInfo struct
			ConsensusMap::FileDescription channel_as_map;
			// label is the channel + description provided in the Params
			channel_as_map.label = "iTRAQ_" + String(cm_it->second.name) + "_" + String(cm_it->second.description);
			channel_as_map.size = ms_exp_MS2.size();
			//TODO what about .filename? leave empty?
			// add some more MetaInfo
			channel_as_map.setMetaValue ("channel_name", cm_it->second.name);
			channel_as_map.setMetaValue ("channel_id", cm_it->second.id);
			channel_as_map.setMetaValue ("channel_description", cm_it->second.description);
			channel_as_map.setMetaValue ("channel_center", cm_it->second.center);
			channel_as_map.setMetaValue ("channel_active", String(cm_it->second.active ? "true" : "false"));
			consensus_map.getFileDescriptions()[index++] = channel_as_map;
		}

		// create consensusElements

		Peak2D::CoordinateType allowed_deviation = (Peak2D::CoordinateType) param_.getValue("reporter_mass_shift");
		// now we have picked data
		// --> assign peaks to channels
		UInt element_index = 0;
		for (MSExperiment<>::ConstIterator it=ms_exp_MS2.begin(); it!=ms_exp_MS2.end(); ++it)
		{
			// store RT&MZ of parent ion as centroid of ConsensusFeature
			ConsensusFeature cf;
			cf.setRT(it->getRT());
			if (it->getPrecursors().size()==1)
			{
				cf.setMZ(it->getPrecursors()[0].getMZ());
			}
			else
			{
				//TODO meckern?
			}

			Peak2D channel_value;
			channel_value.setRT(it->getRT());
			// for each each channel
			Int index = 0;
			Peak2D::IntensityType overall_intensity = 0;
			for (ChannelMapType::const_iterator cm_it = channel_map_.begin(); cm_it!=channel_map_.end(); ++cm_it)
			{
				// set mz-position of channel
				channel_value.setMZ(cm_it->second.center);
				// reset intensity
				channel_value.setIntensity(0);

				//add up all signals
				for (MSExperiment<>::SpectrumType::ConstIterator mz_it = 
								it->MZBegin(cm_it->second.center - allowed_deviation)
						 ;	mz_it != it->MZEnd(cm_it->second.center + allowed_deviation)
						 ;	++mz_it
				)
				{
					channel_value.setIntensity( channel_value.getIntensity() + mz_it->getIntensity());
				}
				
				overall_intensity += channel_value.getIntensity();

				// add channel to ConsensusFeature
				cf.insert (index++, element_index, channel_value);

			} // ! channel_iterator

			
			// check featureHandles are not empty
			if (overall_intensity>0)
			{
				cf.setIntensity(overall_intensity);
				consensus_map.push_back(cf);
			}
			else
			{
				std::cout << "iTRAQ: skipping scan from RT " << cf.getRT() << " & MZ " << cf.getMZ() << "due to no iTRAQ information\n";
			}

			// the tandem-scan in the order they appear in the experiment
			++element_index;
		} // ! Experiment iterator

		#ifdef ITRAQ_DEBUG
		std::cout << "processed " << element_index << " scans" << std::endl;
		#endif

		consensus_map.setExperimentType("itraq");
			
		return;
	}

	void ItraqChannelExtractor::setDefaultParams_()
	{
		defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (left to right) in Da from the expected postion (of e.g. 114.1, 115.1)"); 
		defaults_.setMinFloat ("reporter_mass_shift", 0.00000001);
		defaults_.setMaxFloat ("reporter_mass_shift", 0.5);

		defaults_.setValue("intensity_method", "sum", "which method to use for collecting peaks for each channel", StringList::create("advanced")); 
		defaults_.setValidStrings("intensity_method", StringList::create("sum, peak_picker"));

		if (itraq_type_ == ItraqConstants::FOURPLEX)
		{
			defaults_.setValue("channel_active", StringList::create("114:myReference"), "Each channel that was used in the experiment and its description (114-117 for 4plex; 113-121 for 8-plex) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."); 
		}
		else
		{
			defaults_.setValue("channel_active", StringList::create("113:myReference"), "Each channel that was used in the experiment and its description (114-117 for 4plex; 113-121 for 8-plex) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\"."); 
		}

		PeakPickerCWT ppcwt; 
		defaults_.insert ("PeakPicker:", ppcwt.getDefaults());

		defaultsToParam_();

		return;
	}

	/// implemented for DefaultParamHandler
	void ItraqChannelExtractor::updateMembers_()
	{

		// extract channel names
		StringList channels = StringList(param_.getValue("channel_active"));
		ItraqConstants::updateChannelMap(channels, channel_map_);
	}

	/// initialize
	void ItraqChannelExtractor::init_() 
	{
		ItraqConstants::initChannelMap(itraq_type_, channel_map_);
	}

} // ! namespace
 
