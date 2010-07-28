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
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/CHEMISTRY/AdvancedTheoreticalSpectrumGenerator.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
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
		defaults_.setValue("status", "disabled", "Create Tandem-MS scans?"); 
		defaults_.setValidStrings("status", StringList::create("disabled,precursor,MS^E"));

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
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (ESI or MALDI)");
    defaults_.setValidStrings("ionization_type", StringList::create("ESI,MALDI"));

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

  void RawTandemMSSignalSimulation::generateMSESpectra_(const FeatureMapSim & features, const MSSimExperiment & experiment, MSSimExperiment & ms2)
  {
    AdvancedTheoreticalSpectrumGenerator adv_spec_gen;
    adv_spec_gen.loadProbabilisticModel();
    Param p;
    p.setValue("block_method:rt_block_size", features.size()); // merge all single spectra
    p.setValue("block_method:ms_levels", IntList::create("2"));
    SpectraMerger sm;
    sm.setParameters(p);

    DoubleReal sampling_rate = 1;
    //guess sampling rate from two adjacent full scans:
    if (experiment.size()>=2) sampling_rate = experiment[1].getRT() - experiment[0].getRT();

    // validate features Metavalues exist and are valid:
    for (Size i_f=0;i_f<features.size();++i_f)
    {
      if (!features[i_f].metaValueExists("elution_profile_bounds")
          ||
          !features[i_f].metaValueExists("elution_profile_intensities"))
      {
        throw Exception::ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,"MetaValue:elution_profile_***");
      }
      // check if values fit the experiment:
      const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
      OPENMS_PRECONDITION(elution_bounds[0] < experiment.size(), "Elution profile out of bounds (left)");
      OPENMS_PRECONDITION(elution_bounds[2] < experiment.size(), "Elution profile out of bounds (right)");
      OPENMS_PRECONDITION(experiment[elution_bounds[0]].getRT() == elution_bounds[1], "Elution profile RT shifted (left)");
      OPENMS_PRECONDITION(experiment[elution_bounds[2]].getRT() == elution_bounds[3], "Elution profile RT shifted (right)");
      const DoubleList& elution_ints   = features[i_f].getMetaValue("elution_profile_intensities");
      OPENMS_PRECONDITION(elution_bounds[2] - elution_bounds[0] + 1 == elution_ints.size(), "Elution profile size does not match bounds");
    }

    for (Size i=0;i<experiment.size();++i)
    { // create MS2 for every MS scan
      
      // check which features elute in the current MS scan
      std::vector <Size> features_fragmented;
      for (Size i_f=0;i_f<features.size();++i_f)
      {
        const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
        if ((elution_bounds[1] <= experiment[i].getRT()) && (experiment[i].getRT() <= elution_bounds[3]))
        {
          features_fragmented.push_back(i_f);
        }
      }

      if (features_fragmented.size()==0) continue;

      // now we have all features that elute in this scan -> create MS2 scans
      MSExperiment<RichPeak1D> MS2_spectra;
      MS2_spectra.resize(features_fragmented.size());

      for (Size index=0;index<features_fragmented.size();++index)
      {
        Size i_f = features_fragmented[index];
        // create spectrum
        // todo: we could do this once per feature and not for every scan ... would be a lot faster.. but realistic? Sandro??
        AASequence seq = features[i_f].getPeptideIdentifications()[0].getHits()[0].getSequence();
        adv_spec_gen.simulate(MS2_spectra[index], seq, rnd_gen_,features[i_f].getCharge());
        MS2_spectra[index].setMSLevel(2);
        MS2_spectra[index].setRT(experiment[i].getRT() + sampling_rate*(double(index+1) / double(features_fragmented.size()+2)));
        // adjust intensity of single MS2 spectra by feature intensity
        const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
        const DoubleList& elution_ints   = features[i_f].getMetaValue("elution_profile_intensities");
        DoubleReal factor = elution_ints [i - elution_bounds[0] ];
        for (MSSpectrum<RichPeak1D>::iterator it=MS2_spectra[index].begin();it!=MS2_spectra[index].end();++it)
        {
          it->setIntensity(it->getIntensity() * factor);
        }
      }
      
      // debug: also add single spectra
      for (Size ii=0;ii<MS2_spectra.size();++ii) ms2.push_back(MS2_spectra[ii]); // DEBUG

      // merge all MS2 spectra 
      sm.mergeSpectraBlockWise(MS2_spectra);
      if (MS2_spectra.size()!=1) throw Exception::InvalidSize(__FILE__,__LINE__,__PRETTY_FUNCTION__,MS2_spectra.size() );
      // store merged spectrum
      ms2.push_back(MS2_spectra[0]);

    }

  }

  void RawTandemMSSignalSimulation::generatePrecursorSpectra_(const FeatureMapSim & features, const MSSimExperiment & experiment, MSSimExperiment & ms2)
  {
		IntList qs = (IntList) param_.getValue("Precursor:charge_filter");
		std::set<Int> qs_set(qs.begin(),qs.end());

		//** precursor selection **//
		OfflinePrecursorIonSelection ps;
		Param param = param_.copy("Precursor:",true);
		param.remove("charge_filter");
		ps.setParameters(param);
		// different selection strategies for MALDI and ESI
		bool is_MALDI = (String)param_.getValue("ionization_type") == "MALDI";
    ps.makePrecursorSelectionForKnownLCMSMap(features, experiment,ms2,qs_set,is_MALDI);

		
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
        
  }

  void RawTandemMSSignalSimulation::generateRawTandemSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {
		LOG_INFO << "Tandem MS Simulation ... ";
		
    // will hold the MS2 scans
		MSSimExperiment ms2;

    if (param_.getValue("status") == "disabled")
		{
			LOG_INFO << "disabled\n";
			return;
		}
    else if (param_.getValue("status") == "precursor")
		{
			LOG_INFO << "precursor\n";
      generatePrecursorSpectra_ (features, experiment, ms2);
    }
    else // MS^E
    {
			LOG_INFO << "MS^E\n";
      generateMSESpectra_ (features, experiment, ms2);
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
