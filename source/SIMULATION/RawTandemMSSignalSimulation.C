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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorSet.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/SYSTEM/File.h>


namespace OpenMS
{

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const SimRandomNumberGenerator& rng)
  : DefaultParamHandler("RawTandemMSSignalSimulation"),
    rnd_gen_(&rng)
  {
		// Tandem MS params
		defaults_.setValue("status", "disabled", "Create Tandem-MS scans?"); 
		defaults_.setValidStrings("status", StringList::create("disabled,precursor,MS^E"));

		//TODO: we should think of more ways to select precursors
		defaults_.setValue("Precursor:charge_filter",IntList::create(StringList::create("2,3")), "Charges considered for MS2 fragmentation."); 
		defaults_.setMinInt("Precursor:charge_filter",1);
		defaults_.setMaxInt("Precursor:charge_filter",5);
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
		defaults_.setValue("MS_E:add_single_spectra","false","If true, the MS2 spectra for each peptide signal are included in the output (might be a lot). They will have a meta value 'MS_E_debug_spectra' attached, so they can be filtered out. Full MS_E spectra will have 'MS_E_FullSpectrum' instead.");
		defaults_.setValidStrings("MS_E:add_single_spectra", StringList::create("true,false"));
		defaults_.setValue("tandem_mode", 0, "Algorithm to generate the tandem-MS spectra. 0 - fixed intensities, 1 - SVC prediction (abundant/missing), 2 - SVR prediction of peak intensity \n");
		defaults_.setMinInt("tandem_mode",0);
		defaults_.setMaxInt("tandem_mode",2);
		defaults_.setValue("svm_model_set_file", "examples/simulation/SvmModelSet.dat", "File containing the Filenames of SVM Models for different charge variants");

    subsections_.push_back("TandemSim:");
    defaults_.insert("TandemSim:Simple:", TheoreticalSpectrumGenerator().getDefaults());
    Param svm_par=SvmTheoreticalSpectrumGenerator().getDefaults();
    svm_par.remove("svm_mode");
    svm_par.remove("model_file_name");
    defaults_.insert("TandemSim:SVM:", svm_par);

		// sync'ed Param (also appears in IonizationSimulation)
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaults_.setValidStrings("ionization_type", StringList::create("MALDI,ESI"));

    defaultsToParam_();
  }


  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation& source)
    : DefaultParamHandler(source)
	{
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
  }

  RawTandemMSSignalSimulation& RawTandemMSSignalSimulation::operator = (const RawTandemMSSignalSimulation& source)
  {
    DefaultParamHandler::operator=(source);
		setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    return *this;
  }

  RawTandemMSSignalSimulation::~RawTandemMSSignalSimulation()
  {}

  void RawTandemMSSignalSimulation::generateMSESpectra_(const FeatureMapSim & features, const MSSimExperiment & experiment, MSSimExperiment & ms2)
  {   
    //get tandem mode
    Size tandem_mode=param_.getValue("tandem_mode");

    Param simple_gen_params=param_.copy("TandemSim:Simple:", true);

    TheoreticalSpectrumGenerator simple_generator;
    simple_generator.setParameters(simple_gen_params);

    SvmTheoreticalSpectrumGeneratorSet svm_spec_gen_set;
    //this set will hold the precursor charages that have an Svm model
    std::set<Size>svm_model_charges;

    //if SVR or SVC shall be used
    if(tandem_mode)
    {
      String svm_filename = param_.getValue("svm_model_set_file");
      if (!File::readable( svm_filename ) )
      { // look in OPENMS_DATA_PATH
        svm_filename = File::find( svm_filename );
      }

      svm_spec_gen_set.load(svm_filename);
      svm_spec_gen_set.getSupportedCharges(svm_model_charges);

      //set the parameters for each model
      Param svm_gen_params=param_.copy("TandemSim:SVM:", true);
      std::set<Size>::iterator it;
      for(it=svm_model_charges.begin(); it!=svm_model_charges.end(); ++it)
      {
        svm_spec_gen_set.getSvmModel(*it).setParameters(svm_gen_params);
      }
    }

    Param p;
    p.setValue("block_method:rt_block_size", features.size()); // merge all single spectra
    p.setValue("block_method:ms_levels", IntList::create("2"));
    SpectraMerger sm;
    sm.setParameters(p);

    DoubleReal sampling_rate = 1;
    //guess sampling rate from two adjacent full scans:
    if (experiment.size()>=2) sampling_rate = experiment[1].getRT() - experiment[0].getRT();

    MSSimExperiment single_ms2_spectra;
    single_ms2_spectra.resize(features.size());


    // preparation & validation of input
    for (Size i_f=0;i_f<features.size();++i_f)
    {
      // sample MS2 spectra for each feature
      AASequence seq = features[i_f].getPeptideIdentifications()[0].getHits()[0].getSequence();
      //TODO: work around RichPeak1D restriction
      RichPeakSpectrum tmp_spec;      
      Int prec_charge=features[i_f].getCharge();

      if(tandem_mode && svm_model_charges.count(prec_charge))
      {
        svm_spec_gen_set.simulate(tmp_spec, seq, rnd_gen_->biological_rng, prec_charge);        
      }
      else
      {
        simple_generator.getSpectrum(tmp_spec, seq, prec_charge);        
      }
      for(Size peak=0; peak<tmp_spec.size(); ++peak)
      {
        Peak1D p=tmp_spec[peak];
        single_ms2_spectra[i_f].push_back(p);
      }      

      single_ms2_spectra[i_f].setMSLevel(2);
      Precursor prec;
      prec.setMZ(features[i_f].getMZ());
      single_ms2_spectra[i_f].setPrecursors(std::vector<Precursor>(1,prec));
      single_ms2_spectra[i_f].setMetaValue("FeatureID",(String)features[i_f].getUniqueId());


      // validate features Metavalues exist and are valid:
      if (!features[i_f].metaValueExists("elution_profile_bounds")
          ||
          !features[i_f].metaValueExists("elution_profile_intensities"))
      {
        throw Exception::ElementNotFound(__FILE__,__LINE__,__PRETTY_FUNCTION__,"MetaValue:elution_profile_***");
      }
      // check if values fit the experiment:

      #ifdef OPENMS_ASSERTIONS
      const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
      const DoubleList& elution_ints   = features[i_f].getMetaValue("elution_profile_intensities");

      OPENMS_PRECONDITION(elution_bounds[0] < experiment.size(), "Elution profile out of bounds (left)");
      OPENMS_PRECONDITION(elution_bounds[2] < experiment.size(), "Elution profile out of bounds (right)");
      OPENMS_PRECONDITION(experiment[elution_bounds[0]].getRT() == elution_bounds[1], "Elution profile RT shifted (left)");
      OPENMS_PRECONDITION(experiment[elution_bounds[2]].getRT() == elution_bounds[3], "Elution profile RT shifted (right)");      
      OPENMS_PRECONDITION(elution_bounds[2] - elution_bounds[0] + 1 == elution_ints.size(), "Elution profile size does not match bounds");
      #endif
    }

    // creating the MS^E scan:
    bool add_debug_spectra = ((String)param_.getValue("MS_E:add_single_spectra") == "true");

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
      MSSimExperiment MS2_spectra;
      MS2_spectra.resize(features_fragmented.size());

      for (Size index=0;index<features_fragmented.size();++index)
      {
        Size i_f = features_fragmented[index];
        // create spectrum
        MS2_spectra[index] = single_ms2_spectra[i_f];
        MS2_spectra[index].setRT(experiment[i].getRT() + sampling_rate*(double(index+1) / double(features_fragmented.size()+2)));
        // adjust intensity of single MS2 spectra by feature intensity
        const DoubleList& elution_bounds = features[i_f].getMetaValue("elution_profile_bounds");
        const DoubleList& elution_ints   = features[i_f].getMetaValue("elution_profile_intensities");
        DoubleReal factor = elution_ints [i - elution_bounds[0] ];
        for (MSSimExperiment::SpectrumType::iterator it=MS2_spectra[index].begin();it!=MS2_spectra[index].end();++it)
        {          
          it->setIntensity(it->getIntensity() * factor * features[i_f].getIntensity());          
        }
      }
      
      // debug: also add single spectra
      if (add_debug_spectra)
      {
        for (Size ii=0;ii<MS2_spectra.size();++ii)
        {
          ms2.push_back(MS2_spectra[ii]); // DEBUG
          ms2.back().setMetaValue("MS_E_debug_spectra","true");          
        }
      }

      // merge all MS2 spectra 
      sm.mergeSpectraBlockWise(MS2_spectra);
      if (MS2_spectra.size()!=1) throw Exception::InvalidSize(__FILE__,__LINE__,__PRETTY_FUNCTION__,MS2_spectra.size() );
      // store merged spectrum
      MS2_spectra[0].setMetaValue("MS_E_FullSpectrum","true");
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
		ps.makePrecursorSelectionForKnownLCMSMap(features, experiment, ms2, qs_set, is_MALDI);

		
		//** actual MS2 signal **//
		std::cout << "MS2 features selected: " << ms2.size() << "\n";

    TheoreticalSpectrumGenerator simple_generator;
    Param simple_gen_params=param_.copy("TandemSim:Simple:", true);
    simple_generator.setParameters(simple_gen_params);

    //get tandem mode
    Size tandem_mode=param_.getValue("tandem_mode");

    SvmTheoreticalSpectrumGeneratorSet svm_spec_gen_set;
    //this set will hold the precursor charages that have an Svm model
    std::set<Size>svm_model_charges;

    if(tandem_mode)
    {
      String svm_filename = param_.getValue("svm_model_set_file");
      if (!File::readable( svm_filename ) )
      { // look in OPENMS_DATA_PATH
        svm_filename = File::find( svm_filename );
      }

      svm_spec_gen_set.load(svm_filename);
      svm_spec_gen_set.getSupportedCharges(svm_model_charges);

      //set the parameters for each model
      Param svm_gen_params=param_.copy("TandemSim:SVM:", true);
      std::set<Size>::iterator it;
      for(it=svm_model_charges.begin(); it!=svm_model_charges.end(); ++it)
      {
        svm_spec_gen_set.getSvmModel(*it).setParameters(svm_gen_params);
      }
    }

		for (Size i = 0; i < ms2.size(); ++i)
    {
		  IntList ids = (IntList) ms2[i].getMetaValue("parent_feature_ids");
      DoubleReal prec_intens = ms2[i].getPrecursors()[0].getIntensity();
      MSSimExperiment tmp_spectra;
		  for(Size id =0; id<ids.size();++id)
		  {
		    AASequence seq = features[ids[id]].getPeptideIdentifications()[0].getHits()[0].getSequence();        
        RichPeakSpectrum tmp_spec;
        std::cerr<<"generate Spec for peptide: "<<seq<<std::endl;
        Int prec_charge=features[ids[id]].getCharge();

        if(tandem_mode && svm_model_charges.count(prec_charge))
        {
          svm_spec_gen_set.simulate(tmp_spec, seq, rnd_gen_->biological_rng, prec_charge);
        }
        else
        {
          simple_generator.getSpectrum(tmp_spec, seq, prec_charge);
        }
        //svm_spec_gen.simulate(tmp_spec, seq, rnd_gen_->biological_rng, features[ids[id]].getCharge());
        //TODO: work around RichPeak1D restriction        
        for(Size peak=0; peak<tmp_spec.size(); ++peak)
        {
          Peak1D p=tmp_spec[peak];
          p.setIntensity(p.getIntensity()*prec_intens);
          ms2[i].push_back(p);
        }
      }
      ms2[i].sortByPosition();
    }
  }

  void RawTandemMSSignalSimulation::generateRawTandemSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {
		LOG_INFO << "Tandem MS Simulation ... ";
		
    // will hold the MS2 scans
		MSSimExperiment ms2;

    if (param_.getValue("status") == "disabled")
		{
			LOG_INFO << "disabled" << std::endl;
			return;
		}
    else if (param_.getValue("status") == "precursor")
		{
			LOG_INFO << "precursor" << std::endl;
      generatePrecursorSpectra_ (features, experiment, ms2);
    }
    else // MS^E
    {
			LOG_INFO << "MS^E" << std::endl;
      generateMSESpectra_ (features, experiment, ms2);
    }

		// append MS2 to experiment
		experiment.insert(experiment.end(), ms2.begin(), ms2.end());

  }

  
}
