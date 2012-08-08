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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/KERNEL/RangeUtils.h>

namespace OpenMS
{
  /// Constructor (assuming 4plex)
  ItraqChannelExtractor::ItraqChannelExtractor() :
    DefaultParamHandler("ItraqChannelExtractor"),
    itraq_type_(ItraqConstants::FOURPLEX),
    channel_map_()
  {
    init_();
    setDefaultParams_();
  }

  /// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES)
  ItraqChannelExtractor::ItraqChannelExtractor(Int itraq_type) :
    DefaultParamHandler("ItraqChannelExtractor"),
    itraq_type_(itraq_type),
    channel_map_()
  {
    init_();
    setDefaultParams_();
  }

  /// Constructor with iTRAQ type (from enum ItraqConstants::ITRAQ_TYPES) and param
  ItraqChannelExtractor::ItraqChannelExtractor(Int itraq_type, const Param & param) :
    DefaultParamHandler("ItraqChannelExtractor"),
    itraq_type_(itraq_type),
    channel_map_()
  {
    init_();
    setDefaultParams_();
    setParameters(param);
    updateMembers_();
  }

  /// copy constructor
  ItraqChannelExtractor::ItraqChannelExtractor(const ItraqChannelExtractor & cp) :
    DefaultParamHandler(cp),
    ItraqConstants(cp),
    itraq_type_(cp.itraq_type_),
    channel_map_(cp.channel_map_)
  {

  }

  /// assignment operator
  ItraqChannelExtractor & ItraqChannelExtractor::operator=(const ItraqChannelExtractor & rhs)
  {
    if (this == &rhs)
      return *this;

    DefaultParamHandler::operator=(rhs);
    ItraqConstants::operator=(rhs);
    itraq_type_ = rhs.itraq_type_;
    channel_map_ = rhs.channel_map_;

    return *this;
  }

  /// @brief extracts the iTRAQ channels from the MS data and stores intensity values in a consensus map
  ///
  /// @param ms_exp_data Raw data to read
  /// @param consensus_map Output each MS² scan as a consensus feature
  /// @throws Exception::MissingInformation if no scans present or MS² scan has no precursor
  void ItraqChannelExtractor::run(const MSExperiment<Peak1D> & ms_exp_data, ConsensusMap & consensus_map)
  {
    if (ms_exp_data.empty())
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Experiment has no scans!");
    }

    MSExperiment<> ms_exp_MS2;

    String mode = (String) param_.getValue("select_activation");
    std::cout << "Selecting scans with activation mode: " << (mode == "" ? "any" : mode) << "\n";
    HasActivationMethod<MSExperiment<Peak1D>::SpectrumType> activation_predicate(StringList::create(mode));

    for (size_t idx = 0; idx < ms_exp_data.size(); ++idx)
    {
      if (ms_exp_data[idx].getMSLevel() == 2)
      {
        if (mode == "" || activation_predicate(ms_exp_data[idx]))
        {
          // copy only MS² scans
          ms_exp_MS2.push_back(ms_exp_data[idx]);
        }
        else
        {
          //std::cout << "deleting spectrum # " << idx << " with RT: " << ms_exp_data[idx].getRT() << "\n";
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
    for (ChannelMapType::const_iterator cm_it = channel_map_.begin(); cm_it != channel_map_.end(); ++cm_it)
    {
      // structure of Map cm_it
      //  first == channel-name as Int e.g. 114
      //  second == ChannelInfo struct
      ConsensusMap::FileDescription channel_as_map;
      // label is the channel + description provided in the Params
      if (itraq_type_ != TMT_SIXPLEX)
        channel_as_map.label = "iTRAQ_" + String(cm_it->second.name) + "_" + String(cm_it->second.description);
      else
        channel_as_map.label = "TMT_" + String(cm_it->second.name) + "_" + String(cm_it->second.description);

      channel_as_map.size = ms_exp_MS2.size();
      //TODO what about .filename? leave empty?
      // add some more MetaInfo
      channel_as_map.setMetaValue("channel_name", cm_it->second.name);
      channel_as_map.setMetaValue("channel_id", cm_it->second.id);
      channel_as_map.setMetaValue("channel_description", cm_it->second.description);
      channel_as_map.setMetaValue("channel_center", cm_it->second.center);
      channel_as_map.setMetaValue("channel_active", String(cm_it->second.active ? "true" : "false"));
      consensus_map.getFileDescriptions()[index++] = channel_as_map;
    }

    // create consensusElements

    Peak2D::CoordinateType allowed_deviation = (Peak2D::CoordinateType) param_.getValue("reporter_mass_shift");
    // now we have picked data
    // --> assign peaks to channels
    UInt element_index(0);

    for (MSExperiment<>::ConstIterator it = ms_exp_MS2.begin(); it != ms_exp_MS2.end(); ++it)
    {
      // store RT&MZ of parent ion as centroid of ConsensusFeature
      ConsensusFeature cf;
      cf.setUniqueId();
      cf.setRT(it->getRT());
      if (it->getPrecursors().size() >= 1)
      {
        cf.setMZ(it->getPrecursors()[0].getMZ());
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("No precursor information given for scan native ID ") + String(it->getNativeID()) + " with RT " + String(it->getRT()));
      }

      Peak2D channel_value;
      channel_value.setRT(it->getRT());
      // for each each channel
      Int index = 0;
      Peak2D::IntensityType overall_intensity = 0;
      for (ChannelMapType::const_iterator cm_it = channel_map_.begin(); cm_it != channel_map_.end(); ++cm_it)
      {
        // set mz-position of channel
        channel_value.setMZ(cm_it->second.center);
        // reset intensity
        channel_value.setIntensity(0);

        //add up all signals
        for (MSExperiment<>::SpectrumType::ConstIterator mz_it =
               it->MZBegin(cm_it->second.center - allowed_deviation)
             ; mz_it != it->MZEnd(cm_it->second.center + allowed_deviation)
             ; ++mz_it
             )
        {
          channel_value.setIntensity(channel_value.getIntensity() + mz_it->getIntensity());
        }

        overall_intensity += channel_value.getIntensity();

        // add channel to ConsensusFeature
        cf.insert(index++, channel_value, element_index);

      }       // ! channel_iterator


      // check featureHandles are not empty
      if (overall_intensity == 0)
      {
        cf.setMetaValue("all_empty", String("true"));
      }
      cf.setIntensity(overall_intensity);
      consensus_map.push_back(cf);

      // the tandem-scan in the order they appear in the experiment
      ++element_index;
    }     // ! Experiment iterator


#ifdef ITRAQ_DEBUG
    std::cout << "processed " << element_index << " scans" << std::endl;
#endif

    consensus_map.setExperimentType("itraq");

    return;
  }

  void ItraqChannelExtractor::setDefaultParams_()
  {
    defaults_.setValue("select_activation", Precursor::NamesOfActivationMethod[Precursor::HCID], "Operate only on MSn scans where any of its precursors features a certain activation method (usually HCD for iTRAQ). Set to empty string if you want to disable filtering.\n");
    StringList activation_list(std::vector<std::string>(Precursor::NamesOfActivationMethod, &Precursor::NamesOfActivationMethod[Precursor::SIZE_OF_ACTIVATIONMETHOD - 1]));
    activation_list.push_back("");     // allow disabling this
    defaults_.setValidStrings("select_activation", activation_list);

    defaults_.setValue("reporter_mass_shift", 0.1, "Allowed shift (left to right) in Da from the expected position (of e.g. 114.1, 115.1)");
    defaults_.setMinFloat("reporter_mass_shift", 0.00000001);
    defaults_.setMaxFloat("reporter_mass_shift", 0.5);

    defaults_.setValue("channel_active",
                       (itraq_type_ == TMT_SIXPLEX ? StringList::create("126:liver,131:lung") : StringList::create("114:liver,117:lung")),
                       "Each channel that was used in the experiment and its description (114-117 for 4plex; 113-121 for 8-plex;126-131 for TMT-6-plex) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\".");

    defaultsToParam_();

    return;
  }

  /// implemented for DefaultParamHandler
  void ItraqChannelExtractor::updateMembers_()
  {
    // extract channel names
    ItraqConstants::initChannelMap(itraq_type_, channel_map_);
    ItraqConstants::updateChannelMap(StringList(param_.getValue("channel_active")), channel_map_);
  }

  /// initialize
  void ItraqChannelExtractor::init_()
  {
    ItraqConstants::initChannelMap(itraq_type_, channel_map_);
  }

} // ! namespace
