// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/ID/OfflinePrecursorIonSelection.h>
#include <OpenMS/CHEMISTRY/AdvancedTheoreticalSpectrumGenerator.h>
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
		defaults_.setValue("enabled", "false", "Create Tandem-MS scans?"); 
		defaults_.setValidStrings("enabled", StringList::create("true,false"));

		//TODO: we should think of more ways to select precursors
		defaults_.setValue("Precursor:charge_filter",IntList::create(StringList::create("2,3")), "Charges considered for MS2 fragmentation."); 
		defaults_.setMinInt("Precursor:charge_filter",1);
		defaults_.setMaxInt("Precursor:charge_filter",30);
		defaults_.setValue("Precursor:ms2_spectra_per_rt_bin",5,"Number of allowed MS/MS spectra in a retention time bin.");
		defaults_.setMinInt("Precursor:ms2_spectra_per_rt_bin",1);
		defaults_.setValue("Precursor:exclude_overlapping_peaks","true","If true overlapping or nearby peaks (within min_peak_distance) are excluded for selection.");
		defaults_.setValidStrings("Precursor:exclude_overlapping_peaks", StringList::create("true,false"));
		defaults_.setValue("Precursor:min_peak_distance",2.,"The minimal distance (in Da) of two peaks in one spectrum so that they can be selected.");
		defaults_.setMinFloat("Precursor:min_peak_distance",0.);
		defaults_.setValue("Precursor:use_dynamic_exclusion","true","If true dynamic exclusion is applied.");
		defaults_.setValidStrings("Precursor:use_dynamic_exclusion", StringList::create("true,false"));
		defaults_.setValue("Precursor:exclusion_time",100.,"The time (in seconds) a feature is excluded after its last selection.");
		defaults_.setMinFloat("Precursor:exclusion_time",0.);
		
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
		std::cout << "Tandem MS Simulation ...\n";
		if (param_.getValue("enabled") == "false")
		{
			std::cout << " disabled\n";
			return;
		}

		// will hold the selected precursors
		MSSimExperiment ms2;
		IntList qs = (IntList) param_.getValue("Precursor:charge_filter");
		std::set<Int> qs_set(qs.begin(),qs.end());

		//** precursor selection **//
		OfflinePrecursorIonSelection ps;
		Param param = param_.copy("Precursor:",true);
		param.remove("charge_filter");
		ps.setParameters(param);
		// different selection strategies for MALDI and ESI
		if((String)param_.getValue("ionization_type") == "ESI")
			{
				ps.makePrecursorSelectionForKnownLCMSMap(features, experiment,ms2,qs_set,false);
			}
		else ps.makePrecursorSelectionForKnownLCMSMap(features, experiment,ms2,qs_set,true);
		
		//** actual MS2 signal **//
		std::cout << "MS2 features selected: " << ms2.size() << "\n";
		
		AdvancedTheoreticalSpectrumGenerator adv_spec_gen;
		adv_spec_gen.loadProbabilisticModel();
		for (Size i = 0; i < ms2.size(); ++i)
    {
		  //RichPeakSpectrum ms2_tmp(ms2[i].size());
		  IntList ids = (IntList) ms2[i].getMetaValue("parent_feature_ids");
//		  for(Size pk=0; pk<ms2[i].size();++pk)
//		    ms2_tmp[pk]=ms2[i][pk];
		  for(Size id =0; id<ids.size();++id)
		  {
		    AASequence seq = features[ids[id]].getPeptideIdentifications()[0].getHits()[0].getSequence();
		    adv_spec_gen.simulate(ms2[i], seq, rnd_gen_,1);
      }
//		  for(Size pk=0; pk<ms2[i].size();++pk)
//		    ms2[i][pk]=ms2_tmp[pk];
    }

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
				// reset sum matrix to 0
				gsl_matrix_scale (itraq_intensity_sum, 0);
				
				// add up signal of all features
				// TODO: take care of actual position of feature relative to precursor!
				IntList parent_fs = (IntList) it->getMetaValue("parent_feature_ids");
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
