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
#include <gsl/gsl_blas.h>

namespace OpenMS
{

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const gsl_rng * random_generator)
  : DefaultParamHandler("RawTandemMSSignalSimulation"),
		itraq_type_(),
		channel_map_(),
		isotope_corrections_(),
		rnd_gen_(random_generator)
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
		defaults_.setValue("enabled", "true", "Create Tandem-MS scans?"); 
		defaults_.setValidStrings("enabled", StringList::create("true,false"));

		//TODO: we should think of more ways to select precursors
		defaults_.setValue("Precursor:ChargeFilter",IntList::create(StringList::create("2,3")), "Charges considered for MS2 fragmentation."); 
		defaults_.setMinInt("Precursor:ChargeFilter",1);
		defaults_.setMaxInt("Precursor:ChargeFilter",30);

		// sync'ed Param (also appears in IonizationSimulation)
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", StringList::create("MALDI,ESI"));

		// iTRAQ
		defaults_.setValue("iTRAQ:iTRAQ", "off", "off,4plex or 8plex iTRAQ?");
		defaults_.setValidStrings("iTRAQ:iTRAQ", StringList::create("off,4plex,8plex"));

		defaults_.setValue("iTRAQ:reporter_mass_shift", 0.1, "Allowed shift (uniformly distributed - left to right) in Da from the expected postion (of e.g. 114.1, 115.1)"); 
		defaults_.setMinFloat ("iTRAQ:reporter_mass_shift", 0);
		defaults_.setMaxFloat ("iTRAQ:reporter_mass_shift", 0.5);

		defaults_.setValue("iTRAQ:channel_active_4plex", StringList::create("114:myReference"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."); 
		defaults_.setValue("iTRAQ:channel_active_8plex", StringList::create("113:myReference"), "Eight-plex only: Each channel that was used in the experiment and its description (113-121) in format <channel>:<name>, e.g. \"113:myref\",\"115:liver\",\"118:lung\"."); 

		StringList isotopes = ItraqConstants::getIsotopeMatrixAsStringList(itraq_type_, isotope_corrections_);
		defaults_.setValue("isotope_correction_values", isotopes, "override default values (see Documentation); use the following format: <channel>:<-2Da>/<-1Da>/<+1Da>/<+2Da> ; e.g. '114:0/0.3/4/0' , '116:0.1/0.3/3/0.2' ", StringList::create("advanced"));
		
    defaultsToParam_();
  }

  void RawTandemMSSignalSimulation::updateMembers_()
  {
		StringList channels_active;
		if (param_.getValue("iTRAQ:iTRAQ") == "off") 
		{
			return;
		}
		if (param_.getValue("iTRAQ:iTRAQ") == "4plex") 
		{
			itraq_type_ = ItraqConstants::FOURPLEX;
			channels_active = param_.getValue("iTRAQ:channel_active_4plex");
		}
		else if (param_.getValue("iTRAQ:iTRAQ") == "8plex") 
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
		if (param_.getValue("enabled") == "false") return;

		// will hold the selected precursors
		MSSimExperiment ms2;
		IntList qs = (IntList) param_.getValue("Precursor:ChargeFilter");
		std::set<Int> qs_set(qs.begin(),qs.end());

		//** precursor selection **//

		//TODO: introduce white & blacklists?
		//TODO: introduce some notion of MALDI spot capacity
		
		// current dummy function: naively select 40 highest intensity features
		features.sortByIntensity(true);
		for (Size i=0; i<40 && i<features.size(); ++i)
		{

			// charge not in "ChargeFilter" list
			if (qs_set.count(features[i].getCharge())<1) continue;
			

			MSSimExperiment::iterator scan = experiment.RTBegin(features[i].getRT());
			MSSimExperiment::SpectrumType ms2_spec;
			Precursor p;
			std::vector< Precursor > pcs;
			p.setIntensity(features[i].getIntensity());
			p.setMZ(features[i].getMZ());
			p.setCharge(features[i].getCharge());
			pcs.push_back(p);
			ms2_spec.setPrecursors(pcs);
			ms2_spec.setRT(scan->getRT());
			ms2_spec.setMSLevel(2);
			ms2.push_back(ms2_spec);
			// link ms2 spectrum with features overlapping its precursor
			// Warning: this depends on the current order of features in the map
			// Attention: make sure to name ALL features that overlap, not only one!
			ms2.setMetaValue("parent_feature_ids", IntList::create(String(i)));
			std::cout << " MS2 spectra generated at: " << scan->getRT() << " x " << p.getMZ() << "\n";
		}

		//** actual MS2 signal **//

		// TODO: Sandro, your turn :)


		//** iTRAQ reporters **//
		if (param_.getValue("iTRAQ:iTRAQ") != "off")
		{
		
			std::cout << "Matrix used: \n" << ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_) << "\n\n";
				
			gsl_matrix* channel_frequency = ItraqConstants::translateIsotopeMatrix(itraq_type_, isotope_corrections_).toGslMatrix();
			gsl_matrix* itraq_intensity_observed = Matrix<SimIntensityType>(ItraqConstants::CHANNEL_COUNT[itraq_type_],1).toGslMatrix();
			gsl_matrix* itraq_intensity_sum = Matrix<SimIntensityType>(ItraqConstants::CHANNEL_COUNT[itraq_type_],1).toGslMatrix();
			
			std::vector< Matrix<Int> > channel_names(2);
			channel_names[0].setMatrix<4,1>(ItraqConstants::CHANNELS_FOURPLEX);
			channel_names[1].setMatrix<8,1>(ItraqConstants::CHANNELS_EIGHTPLEX);

			// add signal...
			for (MSSimExperiment::iterator it=ms2.begin(); it!=ms2.end(); ++it)
			{
				std::cout << "adding iTRAQ to MS2 @ " << it->getRT() << "\n";
				
				// reset sum matrix to 0
				gsl_matrix_scale (itraq_intensity_sum, 0);
				
				// add up signal of all features
				// TODO: take care of actual position of feature relative to precursor!
				IntList parent_fs = ms2.getMetaValue("parent_feature_ids");
				for (Size i_f=0; i_f < parent_fs.size(); ++i_f)
				{
					// apply isotope matrix to active channels
					gsl_matrix* row = getItraqIntensity_(features[i_f]).toGslMatrix();
					// row * channel_frequency = observed iTRAQ intensities
					gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
													1.0, channel_frequency, row,
													0.0, itraq_intensity_observed);
          // add result to sum
          gsl_matrix_add (itraq_intensity_sum, itraq_intensity_observed);
					gsl_matrix_free (row);
				}
				
				// add signal to MS2 spectrum
				for (Int i_channel=0; i_channel< ItraqConstants::CHANNEL_COUNT[itraq_type_]; ++i_channel)
				{
					MSSimExperiment::SpectrumType::PeakType p;
					// dummy
					p.setMZ(channel_names[itraq_type_].getValue(i_channel,0) + 0.1);
					p.setIntensity(gsl_matrix_get(itraq_intensity_sum, i_channel, 0));
					std::cout << "inserted iTRAQ peak: " << p << "\n";
					it->push_back(p);
				}
			}
			
			gsl_matrix_free (channel_frequency);
			gsl_matrix_free (itraq_intensity_observed);
			gsl_matrix_free (itraq_intensity_sum);
		}
		

		// append MS2 to experiment
		experiment.insert(experiment.end(), ms2.begin(), ms2.end());

  }
  
  Matrix<SimIntensityType> RawTandemMSSignalSimulation::getItraqIntensity_(const Feature & f) const
  {
		StringList keys;
		f.getKeys(keys);

		// prepare map
		Map <Int, SimIntensityType> channel_intensities;
		std::vector< Matrix<Int> > channel_names(2);
		channel_names[0].setMatrix<4,1>(ItraqConstants::CHANNELS_FOURPLEX);
		channel_names[1].setMatrix<8,1>(ItraqConstants::CHANNELS_EIGHTPLEX);
		for (Int i=0; i<ItraqConstants::CHANNEL_COUNT[itraq_type_]; ++i)
		{
			channel_intensities[channel_names[itraq_type_].getValue(i,0)] = 0;
		}		
		
		// fill map with values present (all missing ones remain 0)
		for (StringList::const_iterator it_key = keys.begin(); it_key != keys.end(); it_key++)
		{
			if (!it_key->hasPrefix("intensity_itraq")) continue;
			Int ch = it_key->substr(String("intensity_itraq").length()).toInt();
			channel_intensities[ch] = f.getMetaValue(*it_key);
			std::cout << "raw itraq intensity: " << ch << "->" << f.getMetaValue(*it_key) << "\n";
		}	
		
		// fill the matrix
		Matrix<SimIntensityType> m(ItraqConstants::CHANNEL_COUNT[itraq_type_], 1, 0);
		Size index=0;
		for (Map <Int, SimIntensityType>::const_iterator it=channel_intensities.begin(); it!=channel_intensities.end(); ++it)
		{
			m.setValue(index ,0, it->second);
			++index;
		}

		return m;
		
  }
  
}
