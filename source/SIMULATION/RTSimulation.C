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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RTSimulation.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;


namespace OpenMS {

  
  RTSimulation::RTSimulation(const SimRandomNumberGenerator& random_generator)
    : DefaultParamHandler("RTSimulation"), rnd_gen_(&random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }
  
  
  RTSimulation::RTSimulation(const RTSimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RTSimulation& RTSimulation::operator = (const RTSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }
  
  RTSimulation::~RTSimulation()
  {
  }
  
  void RTSimulation::setDefaultParams_() 
  {
		defaults_.setValue("rt_column", "HPLC", "Modelling of an RT or CE column");
    defaults_.setValidStrings("rt_column", StringList::create("none,HPLC,CE"));
    
    // scaling
    defaults_.setValue("auto_scale","true","Scale predicted RT's/MT's to given 'total_gradient_time'? If 'true', for CE this means that 'CE:lenght_d', 'CE:length_total', 'CE:voltage' have no influence.");
    defaults_.setValidStrings("auto_scale", StringList::create("true,false"));

		// column settings
    defaults_.setValue("total_gradient_time",2500.0,"The duration [s] of the gradient.");
    defaults_.setMinFloat("total_gradient_time", 0.00001);

    // rt scan window
    defaults_.setValue("scan_window:min",500.0,"Start of RT Scan Window [s]");
    defaults_.setMinFloat("scan_window:min",0);
    defaults_.setValue("scan_window:max",1500.0,"End of RT Scan Window [s]");
    defaults_.setMinFloat("scan_window:max",1);

    // rt spacing
    defaults_.setValue("sampling_rate", 2.0, "Time interval [s] between consecutive scans");
    defaults_.setMinFloat("sampling_rate",0.01);
    defaults_.setMaxFloat("sampling_rate",60.0);

    // rt error
    defaults_.setValue("variation:feature_stddev",3,"Standard deviation of shift in retention time [s] from predicted model (applied to every single feature independently)");
    defaults_.setValue("variation:affine_offset",0,"Global offset in retention time [s] from predicted model");
    defaults_.setValue("variation:affine_scale",1,"Global scaling in retention time from predicted model");
    defaults_.setSectionDescription("variation","Random component that simulates technical/biological variation");

    defaults_.setValue("column_condition:distortion", 0, "Distortion of the elution profiles. Good presets are 0 for a perfect elution profile, 1 for a slightly distorted elution profile etc... For trapping instruments (e.g. Orbitrap) distortion should be >4.");
    defaults_.setMinInt("column_condition:distortion", 0);
    defaults_.setMaxInt("column_condition:distortion", 10);


    defaults_.setValue("profile_shape:width:value", 9, "Width of the Exponential Gaussian Hybrid distribution shape of the elution profile. This does not correspond directly to the width in [s].");
    defaults_.setValue("profile_shape:width:variance", 1.6, "Random component of the width (set to 0 to disable randomness), i.e. scale parameter for the lorentzian variation of the variance (Note: The scale parameter has to be >= 0).");
    defaults_.setMinFloat("profile_shape:width:variance",0.0);
    defaults_.setSectionDescription("profile_shape:width","Width of the EGH elution shape, i.e. the sigma^2 parameter, which is computed using 'value' + rnd_cauchy('variance')");

    defaults_.setValue("profile_shape:skewness:value", 0.1, "Asymmetric component of the EGH. Higher absolute(!) values lead to more skewness (negative values cause fronting, positive values cause tailing). Tau parameter of the EGH, i.e. time constant of the exponential decay of the Exponential Gaussian Hybrid distribution shape of the elution profile.");
    defaults_.setValue("profile_shape:skewness:variance", 0.3, "Random component of skewness (set to 0 to disable randomness), i.e. scale parameter for the lorentzian variation of the time constant (Note: The scale parameter has to be > 0).");
    defaults_.setMinFloat("profile_shape:skewness:variance",0.0);
    defaults_.setSectionDescription("profile_shape:skewness","Skewness of the EGH elution shape, i.e. the tau parameter, which is computed using 'value' + rnd_cauchy('variance')");

    // HPLC specific Parameters
    defaults_.setValue("HPLC:model_file","examples/simulation/RTPredict.model","SVM model for retention time prediction");

    // CE specific Parameters
    defaults_.setValue("CE:pH",3.0,"pH of buffer");
    defaults_.setMinFloat("CE:pH", 0);
    defaults_.setMaxFloat("CE:pH", 14);
    
    defaults_.setValue("CE:alpha",0.5,"Exponent Alpha used to calculate mobility");
    defaults_.setMinFloat("CE:alpha", 0);
    defaults_.setMaxFloat("CE:alpha", 1);

    defaults_.setValue("CE:mu_eo",0.0,"Electroosmotic flow");
    defaults_.setMinFloat("CE:mu_eo", 0);
    defaults_.setMaxFloat("CE:mu_eo", 5);
        
    defaults_.setValue("CE:lenght_d", 70.0 ,"Length of capillary [cm] from injection site to MS");
    defaults_.setMinFloat("CE:lenght_d", 0);
    defaults_.setMaxFloat("CE:lenght_d", 1000);

    defaults_.setValue("CE:length_total", 75.0 ,"Total length of capillary [cm]");
    defaults_.setMinFloat("CE:length_total", 0);
    defaults_.setMaxFloat("CE:length_total", 1000);
    
    defaults_.setValue("CE:voltage", 1000.0 ,"Voltage applied to capillary");
    defaults_.setMinFloat("CE:voltage", 0);

		defaultsToParam_();
  }


  void RTSimulation::updateMembers_()
  {
		rt_model_file_ = param_.getValue("HPLC:model_file");
		if (! File::readable( rt_model_file_ ) )
    { // look in OPENMS_DATA_PATH
      rt_model_file_ = File::find( rt_model_file_ );
    }
		total_gradient_time_ = param_.getValue("total_gradient_time");
		gradient_min_ = param_.getValue("scan_window:min");
		gradient_max_ = param_.getValue("scan_window:max");
    if(gradient_max_ > total_gradient_time_)
    {
      LOG_WARN << "total_gradient_time_ smaller than scan_window:max -> invalid parameters!" << endl;
    }
    
    rt_sampling_rate_    = param_.getValue("sampling_rate");

    egh_variance_location_  = param_.getValue("profile_shape:width:value");
    egh_variance_scale_     = param_.getValue("profile_shape:width:variance");
    if (egh_variance_scale_ < 0.0)
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"The scale parameter for the lorentzian variation of the variance has to be >= 0.");
    }

    egh_tau_location_    = param_.getValue("profile_shape:skewness:value");
    egh_tau_scale_       = param_.getValue("profile_shape:skewness:variance");
    if (egh_tau_scale_ < 0.0)
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"The scale parameter for the lorentzian variation of the time constant has to be >= 0.");
    }

  }
  
  void RTSimulation::noRTColumn_(FeatureMapSim & features)
  {
    for(FeatureMapSim::iterator it_f = features.begin(); it_f != features.end();
        ++it_f)
    {
      (*it_f).setRT(-1);
    }
  }


  /**
   @brief Gets a feature map containing the peptides and predicts for those the retention times
   */
  void RTSimulation::predictRT(FeatureMapSim & features)
  {
    LOG_INFO << "RT Simulation ... started" << std::endl;

    vector< DoubleReal>  predicted_retention_times;
		bool is_relative = (param_.getValue("auto_scale")=="true");
    if (param_.getValue("rt_column") == "none")
    {
			noRTColumn_(features);
			return;
		}
		// CE or HPLC:
		else if (param_.getValue("rt_column") == "CE")
		{
			calculateMT_(features, predicted_retention_times);
		}
		else if (param_.getValue("rt_column") == "HPLC")
		{
			vector< AASequence > peptides_aa_vector(features.size());
			for (Size i = 0; i < features.size(); ++i)
			{
				peptides_aa_vector[i] = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence();
			}
			wrapSVM(peptides_aa_vector, predicted_retention_times);
		}	
    
    // rt error dicing
    SimCoordinateType rt_offset = param_.getValue("variation:affine_offset");
    SimCoordinateType rt_scale  = param_.getValue("variation:affine_scale");
    SimCoordinateType rt_ft_stddev = param_.getValue("variation:feature_stddev");      
    
    FeatureMapSim fm_tmp (features);
    fm_tmp.clear(false);
    StringList deleted_features;
    for (Size i = 0; i < predicted_retention_times.size(); ++i)
    {
      // relative -> absolute RT's (with border)
      if (is_relative)
      {
	      predicted_retention_times[i] *= total_gradient_time_;
	    }
      
      //overwrite RT (if given by user)
      if (features[i].metaValueExists("rt"))
      {
        predicted_retention_times[i] = features[i].getMetaValue("rt");
      }
      // add variation
      SimCoordinateType rt_error = gsl_ran_gaussian(rnd_gen_->technical_rng, rt_ft_stddev) + rt_offset;
      predicted_retention_times[i] = predicted_retention_times[i]*rt_scale + rt_error;
      //overwrite RT [no randomization] (if given by user)
      if (features[i].metaValueExists("RT"))
      {
        predicted_retention_times[i] = features[i].getMetaValue("RT");
      }

      // remove invalid peptides & (later) display removed ones
      if (
          (predicted_retention_times[i] < 0.0) || // check for invalid RT
          (predicted_retention_times[i] > gradient_max_) ||  // check if RT is not in scan window
          (predicted_retention_times[i] < gradient_min_)     // check if RT is not in scan window
          )
      {
        deleted_features.push_back(features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString() + " [" +
                                   String::number(predicted_retention_times[i],2)
                                   + "]");
				continue;
      }
      
      features[i].setRT(predicted_retention_times[i]);

      // determine shape parameters for EGH
      DoubleReal variance = egh_variance_location_ + (egh_variance_scale_==0 ? 0 : gsl_ran_cauchy(rnd_gen_->technical_rng, egh_variance_scale_));
      DoubleReal tau = egh_tau_location_ + (egh_tau_scale_==0? 0 : gsl_ran_cauchy(rnd_gen_->technical_rng, egh_tau_scale_));;

      // resample variance if it is below 0
      // try this only 10 times to avoid endless loop in case of
      // a bad parameter combination
      Size retry_variance_sampling=0;
      while(variance <= 0 && retry_variance_sampling < 9)
      {
        variance = egh_variance_location_ + (egh_variance_scale_==0 ? 0 : gsl_ran_cauchy(rnd_gen_->technical_rng, egh_variance_scale_));
        ++retry_variance_sampling;
      }

      if (variance<=0)
      {
        LOG_ERROR << "Sigma^2 was negative, resulting in a feature with width=0. Tried to resample 10 times and then stopped. Skipping feature!" << std::endl;
        deleted_features.push_back(features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString() + " [" +
                                   String::number(predicted_retention_times[i],2)
                                   + "]");
        continue;
      }

      features[i].setMetaValue("RT_egh_variance", variance);
      features[i].setMetaValue("RT_egh_tau", tau);

      fm_tmp.push_back(features[i]);
    }
    
    // print invalid features:
    if(deleted_features.size() > 0)
    {
      LOG_WARN << "RT prediction gave 'invalid' results for " << deleted_features.size() << " peptide(s), making them unobservable.\n";
      if (deleted_features.size()<100) LOG_WARN << "  " << deleted_features.concatenate("\n  ") << std::endl;
      else LOG_WARN << "  (List is too big to show)" << std::endl;
    }
    // only retain valid features:	
    features.swap(fm_tmp);
    
    features.sortByPosition();
    features.updateRanges ();
    
  }


	/// PKA values as given in Rickard1991
	void RTSimulation::getChargeContribution_(Map< String, double> & q_cterm, 
																	 Map< String, double> & q_nterm,
																	 Map< String, double> & q_aa_basic,
																	 Map< String, double> & q_aa_acidic)
	{
		// the actual constants from the paper:
		String aas = "ARNDCQEGHILKMFPSTWYVBZ";
		const double cterm_pkas[] = { 3.20, 3.20, 2.75, 2.75, 2.75, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 3.20, 2.75, 3.20 };
		const double nterm_pkas[] = { 8.20, 8.20, 7.30, 8.60, 7.30, 7.70, 8.20,8.20,8.20,8.20,8.20, 7.70, 9.20, 7.7, 9.00, 7.30, 8.20, 8.20, 7.70, 8.20, 8.03, 8.00};

		String aa_basic = "HRK";
		double aa_basic_pkas[] = {6.20, 12.50, 10.30};

		String aa_acidic = "DECY";
		double aa_acidic_pkas[] = {3.50, 4.50, 10.30, 10.30};

		// clear target structures
		q_cterm.clear();
		q_nterm.clear();
		q_aa_basic.clear();
		q_aa_acidic.clear();
		
		// get params
		DoubleReal ph = param_.getValue("CE:pH");

		// calculate charges according to constants and conditions:
		
		// C&N term
		for (Size i = 0; i<aas.size(); ++i)
		{
			q_nterm[aas[i]] = + 1/( 1+std::pow(10, +ph - nterm_pkas[i] ));
			q_cterm[aas[i]] = - 1/( 1+std::pow(10, -ph + cterm_pkas[i] ));
		}
		
		// basic AA's
		for (Size i = 0; i<aa_basic.size(); ++i)
		{
			q_aa_basic[aa_basic[i]]   = + 1/( 1+std::pow(10, +ph - aa_basic_pkas[i] ));
		}

		// acidic AA's
		for (Size i = 0; i<aa_acidic.size(); ++i)
		{
			q_aa_acidic[aa_acidic[i]] = - 1/( 1+std::pow(10, -ph + aa_acidic_pkas[i] ));
		}
		// add values for ambigous AA according to dayhoff frequencies
		q_aa_acidic["B"] = q_aa_acidic["D"]*(5.5/(5.5+4.3)) + 0*(4.3/(5.5+4.3)); // D~5.5; N~4.3
		q_aa_acidic["Z"] = q_aa_acidic["E"]*(6.0/(6.0+3.9)) + 0*(3.9/(6.0+3.9)); // E~6.0; Q~3.9
	}


  void RTSimulation::calculateMT_(FeatureMapSim & features, std::vector<DoubleReal>& predicted_retention_times)
  {
		Map< String, double> q_cterm, q_nterm, q_aa_basic, q_aa_acidic;
		getChargeContribution_(q_cterm, q_nterm, q_aa_basic, q_aa_acidic);
		
		DoubleReal alpha = param_.getValue("CE:alpha");
		bool auto_scale = (param_.getValue("auto_scale")=="true");
		DoubleReal c = (auto_scale ? 1 : (DoubleReal)param_.getValue("CE:lenght_d") * (DoubleReal)param_.getValue("CE:length_total") / (DoubleReal)param_.getValue("CE:voltage"));

		predicted_retention_times.resize(features.size());
		
		for (Size i = 0; i < features.size(); ++i)
		{
			String seq = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
			
			// ** determine charge of peptide ** 
			
			DoubleReal charge = 0;	
			// C&N term charge contribution
			if (q_nterm.has(seq[0])) charge +=  q_nterm[seq[0]];
			if (q_cterm.has(seq.suffix(1))) charge +=  q_cterm[seq.suffix(1)];
			
			// sidechains ...
			Map< String, Size > frequency_table;
			features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().getAAFrequencies (frequency_table);
			for (Map< String, Size >::const_iterator it = frequency_table.begin(); it != frequency_table.end(); ++it)
			{
				if (q_aa_basic.has(it->first)) charge +=  q_aa_basic[it->first] * it->second;
				if (q_aa_acidic.has(it->first)) charge +=  q_aa_acidic[it->first] * it->second;
			}
	
			// ** determine mass of peptide
			DoubleReal mass = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().getFormula().getAverageWeight();
	
			// ** mobility (mu = mu_ep + mu_eo = (q/MW^alpha) + mu_eo
			DoubleReal mu = ( charge / std::pow(mass, alpha) ) + (auto_scale ? 0 : (DoubleReal)param_.getValue("CE:mu_eo"));
			predicted_retention_times[i] = c / mu; // this is L_d*L_t / (mu * V)  as "c = L_d*L_t/V"
    }
		
		// ** only when Auto-Scaling is active ** /
		std::vector<DoubleReal> rt_sorted(predicted_retention_times);
		std::sort(rt_sorted.begin(), rt_sorted.end() );

    DoubleReal max_rt = rt_sorted.back();
    
    if (auto_scale)
		{
			max_rt = 1; // highest will be scaled to 1

      //std::cerr << "minRT: " << rt_sorted[0] << "   max: " << rt_sorted.back() << "\n";
			// normalize to 5th - 95th percentile (we want to avoid that few outliers with huge/small MT can compress the others to a small MT range):
      DoubleReal mt_5p = rt_sorted[rt_sorted.size()*5/100];
      DoubleReal mt_95p = rt_sorted[rt_sorted.size()*95/100];
			// ... assume 95% MT range at 95th percentile
			DoubleReal range = std::max(1.0, (mt_95p - mt_5p)*0.9);
			
      //std::cerr << " 5% MT: " << mt_5p << ",   95% MT: " << mt_95p << " Range: " << range << "\n";

      DoubleReal new_offset = mt_5p - range*0.05;

			// scale MT's between 0 and 1 (except for outliers --> which will get <0 or >1)
			for (Size i = 0; i < features.size(); ++i)
			{
				predicted_retention_times[i] = (predicted_retention_times[i] - new_offset) / range;
			}
		}
    
    // the width factor is 1.0 at MT=0 and reaches its max (default 2.0) at MT=max
    DoubleReal rt_widening_max = 2.0;
    for (Size i = 0; i < features.size(); ++i)
		{
      features[i].setMetaValue("RT_CE_width_factor", (predicted_retention_times[i]/max_rt * (rt_widening_max-1) + 1));
    }
				
  }

	void RTSimulation::wrapSVM(std::vector<AASequence>& peptide_sequences,std::vector<DoubleReal>& predicted_retention_times)
	{
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    SVMWrapper svm;
		LibSVMEncoder encoder;
		svm_problem* training_data = NULL;
		SVMData prediction_samples;
		SVMData training_samples;
    UInt k_mer_length = 0;
    DoubleReal sigma = 0.0;
    UInt border_length = 0;
    Size max_number_of_peptides(2000); // hard coding pediction bins; larger values only take more memory, result is not affected

		LOG_INFO << "Predicting RT ... ";
    
    svm.loadModel(rt_model_file_);
    
    // load additional parameters
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      String add_paramfile = rt_model_file_ + "_additional_parameters";
      if (! File::readable( add_paramfile ) )
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation: SVM parameter file " + add_paramfile + " is not readable");
      }
      
      Param additional_parameters;
      additional_parameters.load(add_paramfile);
      
      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
          && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation: No border length defined in additional parameters file.");
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
          && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation: No k-mer length defined in additional parameters file.");
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
      
      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
          && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation: No sigma defined in additional parameters file.");
      }
      
      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
    }
    
    svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
    svm.setParameter(SVMWrapper::SIGMA, sigma);
		
    // loading model data
    String sample_file = rt_model_file_ + "_samples";
    if (! File::readable( sample_file ) )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation: SVM sample file " + sample_file + " is not readable");
    }
    training_samples.load(sample_file);
		svm.setTrainingSample(training_samples);
    svm.setTrainingSample(training_data);

		// use maximally max_number_of_peptides peptide sequence at once
		Size tmp_count = 0;
		Size count = 0;
		std::vector<AASequence>::iterator pep_iter_start = peptide_sequences.begin();
		std::vector<AASequence>::iterator pep_iter_stop = peptide_sequences.begin();
		while(count < peptide_sequences.size())
		{
			while(pep_iter_stop != peptide_sequences.end() && tmp_count < max_number_of_peptides)
			{
				++tmp_count;
				++pep_iter_stop;
			}
			std::vector<AASequence> tmp_peptide_seqs;
			tmp_peptide_seqs.insert(tmp_peptide_seqs.end(),pep_iter_start,pep_iter_stop);
			std::vector<DoubleReal> tmp_rts(tmp_peptide_seqs.size(), 0);
			std::vector<DoubleReal> tmp_pred_rts;
			// Encoding test data
			encoder.encodeProblemWithOligoBorderVectors(tmp_peptide_seqs,k_mer_length, allowed_amino_acid_characters, border_length,prediction_samples.sequences);
			prediction_samples.labels = tmp_rts;
		
			svm.predict(prediction_samples, tmp_pred_rts);
			predicted_retention_times.insert(predicted_retention_times.end(),tmp_pred_rts.begin(),tmp_pred_rts.end());
			pep_iter_start = pep_iter_stop;				
			count += tmp_count;
			tmp_count = 0;
		}
    LibSVMEncoder::destroyProblem(training_data);
		
    LOG_INFO << "done" << endl;
	}

  void RTSimulation::predictContaminantsRT(FeatureMapSim & contaminants)
  {
    // iterate of feature map
    for (Size i = 0; i < contaminants.size(); ++i)
    {

      // assign random retention time
      SimCoordinateType retention_time = gsl_ran_flat(rnd_gen_->technical_rng, 0, total_gradient_time_);
      contaminants[i].setRT(retention_time);
    }
  }
  
  
  bool RTSimulation::isRTColumnOn() const
  {
    return (param_.getValue("rt_column") != "none");
  }
  
  SimCoordinateType RTSimulation::getGradientTime() const
  {
    return total_gradient_time_;
  }
  
  void RTSimulation::createExperiment(MSSimExperiment & experiment)
  {
    // this is a closed intervall (it includes gradient_min_ and gradient_max_)
    Size number_of_scans = Size( (gradient_max_ - gradient_min_) / rt_sampling_rate_) + 1;

    experiment = MSSimExperiment();

    if (isRTColumnOn())
    {
      LOG_INFO << "Creating experiment with #" << number_of_scans << " scans ... ";

      experiment.resize(number_of_scans);

      DoubleReal current_scan_rt = gradient_min_;
      Size id = 1;
      for(MSSimExperiment::iterator exp_it = experiment.begin();
          exp_it != experiment.end();
          ++exp_it)
      {
        (*exp_it).setRT(current_scan_rt);

        String spec_id = String("spectrum=") + id;
        ++id;
        (*exp_it).setNativeID(spec_id);

        // dice & store distortion
        (*exp_it).setMetaValue("distortion", 1);

        // TODO (for CE) store peak broadening parameter
        current_scan_rt += rt_sampling_rate_;
      }

      // smooth the distortion with a moving average filter of width 3.0
      smoothRTDistortion_(experiment);
    }
    else
    {
      LOG_INFO << "Creating experiment with a single scan ... ";

      experiment.resize(1);
      experiment[0].setRT(-1);
      experiment[0].setNativeID("spectrum=1");
    }
    experiment.updateRanges();
    LOG_INFO << "done\n";
  }

//#define MSSIM_DEBUG_MOV_AVG_FILTER

  void RTSimulation::smoothRTDistortion_(MSSimExperiment & experiment)
  {
    // how often do we move over the distortions
    const UInt filter_iterations = param_.getValue("column_condition:distortion");

    DoubleReal previous,current,next;

    for(UInt fi = 0 ; fi < filter_iterations ; ++fi )
    {
      // initialize the previous value on position 0
      previous = (DoubleReal) experiment[0].getMetaValue("distortion");
#ifdef MSSIM_DEBUG_MOV_AVG_FILTER
      LOG_WARN << "d <- c(" << previous << ", ";
      vector< DoubleReal > tmp;
#endif
      for(Size scan = 1 ; scan < experiment.size() - 1 ; ++scan)
      {
        current = (DoubleReal) experiment[scan].getMetaValue("distortion");
        next = (DoubleReal) experiment[scan + 1].getMetaValue("distortion");

        DoubleReal smoothed = (previous + current + next) / 3.0;
        smoothed *= gsl_ran_flat (rnd_gen_->technical_rng, 1.0-std::pow(fi+1.0,2)*0.01, 1.0+std::pow(fi+1.0,2)*0.01); // distortion gets worse round by round
        previous = current;

#ifdef MSSIM_DEBUG_MOV_AVG_FILTER
        LOG_WARN << current << ", ";
        tmp.push_back(smoothed);
#endif
        experiment[scan].setMetaValue("distortion", smoothed);
      }

#ifdef MSSIM_DEBUG_MOV_AVG_FILTER
      LOG_WARN << next << ");" << endl;
      LOG_WARN << "smoothed <- c(";
      LOG_WARN << (DoubleReal) experiment[0].getMetaValue("distortion") << ", ";
      for(Size i = 0 ; i  < tmp.size() ; ++i) {
        LOG_WARN << tmp[i] << ", ";
      }
      LOG_WARN << next << ");" << endl;
#endif
    }
  }

}
