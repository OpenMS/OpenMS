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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>

namespace OpenMS {

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const gsl_rng * random_generator)
  : DefaultParamHandler("RawTandemMSSignalSimulation"), rnd_gen_(random_generator),
		itraq_type_(),
		channel_map_(),
		isotope_corrections_()
  {
    init_();
  }


  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation& source)
    : DefaultParamHandler(source),
			itraq_type_(source.itraq_type_),
			channel_map_(source.channel_map_),
			isotope_corrections_(source.isotope_corrections_)
	{
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RawTandemMSSignalSimulation& RawTandemMSSignalSimulation::operator = (const RawTandemMSSignalSimulation& source)
  {
    DefaultParamHandler::operator=(source);
		itraq_type_ = source.itraq_type_;
		channel_map_ = source.channel_map_;
		isotope_corrections_ = source.isotope_corrections_;
		setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }

  RawTandemMSSignalSimulation::~RawTandemMSSignalSimulation()
  {}

	void RawTandemMSSignalSimulation::init_()
	{
		// this needs to come first!
		isotope_corrections_.resize(2);
		isotope_corrections_[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
		isotope_corrections_[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);

		setDefaultParams_();
    updateMembers_();

	}

  void RawTandemMSSignalSimulation::setDefaultParams_()
  {

		// Tandem MS params

		//TODO: we should think of more ways to select precursors
		defaults_.setValue("PrecursorSelectionCount", 4, "Maximal number of precursors selected in each scan."); 
		defaults_.setMinFloat ("PrecursorSelectionCount", 1);
		defaults_.setMaxFloat ("PrecursorSelectionCount", 10);


		// iTRAQ
		defaults_.setValue("iTRAQ:enable_iTRAQ", "false", "enable iTRAQ reporter ion simulation");
		defaults_.setValidStrings("iTRAQ:enable_iTRAQ", StringList::create("true,false"));

		defaults_.setValue("iTRAQ:reporter_mass_shift", 0.1, "Allowed shift (uniformly distributed - left to right) in Da from the expected postion (of e.g. 114.1, 115.1)"); 
		defaults_.setMinFloat ("iTRAQ:reporter_mass_shift", 0);
		defaults_.setMaxFloat ("iTRAQ:reporter_mass_shift", 0.5);

		defaults_.setValue("iTRAQ:iTRAQ_type", "4plex", "4plex or 8plex iTRAQ?");
		defaults_.setValidStrings("iTRAQ:iTRAQ_type", StringList::create("4plex,8plex"));
		defaults_.setValue("iTRAQ:channel_active_4plex", StringList::create("114:myReference"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."); 
		defaults_.setValue("iTRAQ:channel_active_8plex", StringList::create("113:myReference"), "Eight-plex only: Each channel that was used in the experiment and its description (113-121) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\"."); 

		StringList isotopes = ItraqConstants::getIsotopeMatrixAsStringList(itraq_type_, isotope_corrections_);
		defaults_.setValue("isotope_correction_values", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", StringList::create("advanced"));
		
    defaultsToParam_();
  }

  void RawTandemMSSignalSimulation::updateMembers_()
  {
		StringList channels_active;
		if (param_.getValue("iTRAQ:iTRAQ_type") == "4plex") 
		{
			itraq_type_ = ItraqConstants::FOURPLEX;
			channels_active = param_.getValue("iTRAQ:channel_active_4plex");
		}
		else
		{
			itraq_type_ = ItraqConstants::EIGHTPLEX;
			channels_active = param_.getValue("iTRAQ:channel_active_8plex");
		}

		ItraqConstants::initChannelMap(itraq_type_, channel_map_);
		ItraqConstants::updateChannelMap(channels_active, channel_map_);


		// update isotope_corrections_ Matrix with custom values
		StringList channels = param_.getValue("isotope_correction_values");
		if (channels.size()>0)
		{
			ItraqConstants::updateIsotopeMatrixFromStringList(itraq_type_, channels, isotope_corrections_);
		}

	}

  void RawTandemMSSignalSimulation::generateRawTandemSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {

		// precursor selection:


		// TODO: Sandro, your turn :)


		// iTRAQ stuff
		if (param_.getValue("iTRAQ:enable_iTRAQ") == "true")
		{
				// apply isotope matrix to active channels

				// add signal...

		}
		
  }
  
}
