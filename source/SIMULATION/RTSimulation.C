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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RTSimulation.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace OpenMS {
  /// class constants
  const DoubleReal RTSimulation::gradient_front_offset_ = 400.0;
  const DoubleReal RTSimulation::gradient_total_offset_ = 800.0;  
  
  RTSimulation::RTSimulation(const gsl_rng * random_generator)
    : DefaultParamHandler("RTSimulation"), rnd_gen_(random_generator)
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
  {}
  
  /**
   @brief Gets a feature map containing the peptides and predicts for those the retention times
   */
  void RTSimulation::predictRT(FeatureMapSim & features, MSSimExperiment & experiment)
  {
    createExperiment_(experiment);

    if (isRTColumnOn())
    {
      predictFeatureRT_(features);
    }
    else
    {
      noRTColumn_(features);
    }

    for(FeatureMapSim::iterator it_f = features.begin(); it_f != features.end();
        ++it_f)
    {
			double symmetry = gsl_ran_flat (rnd_gen_, symmetry_down_, symmetry_up_);
			double width = gsl_ran_flat (rnd_gen_, 5, 15);
      it_f->setMetaValue("rt_symmetry", symmetry);
      it_f->setMetaValue("rt_width", width);
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
  void RTSimulation::predictFeatureRT_(FeatureMapSim & features)
  {
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    SVMWrapper svm;
    LibSVMEncoder encoder;
    vector<DoubleReal> predicted_retention_times;
    
    svm_problem* training_data = NULL;
    svm_problem* prediction_data = NULL;
    
    UInt k_mer_length = 0;
    DoubleReal sigma = 0.0;
    UInt border_length = 0;
    
		std::cout << "Predicting RT..    " << endl;
    
    // not that elegant...
    vector< String > peptidesVector(features.size());

    for (Size i = 0; i < features.size(); ++i)
    {
      peptidesVector[i] = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
    }
    
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
    
    // Encoding test data
    vector<DoubleReal> rts;
    rts.resize(peptidesVector.size(), 0);
    prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptidesVector, rts, k_mer_length, allowed_amino_acid_characters, border_length);
    
    // loading training data
    String sample_file = rt_model_file_ + "_samples";
    if (! File::readable( sample_file ) )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation: SVM sample file " + sample_file + " is not readable");
    }
    training_data = encoder.loadLibSVMProblem(sample_file);
    
    svm.setTrainingSample(training_data);
    svm.predict(prediction_data, predicted_retention_times);
    
    cout << "Done." << endl;
    delete training_data;
    delete prediction_data;
     
    /// rt error stuff
    SimCoordinateType rt_shift_mean  = param_.getValue("rt_shift_mean");
    SimCoordinateType rt_shift_stddev = param_.getValue("rt_shift_stddev");      
    
    for (Size i = 0; i < peptidesVector.size(); ++i)
    {
      // TODO: remove those peptides + debug out removed ones??
      if (predicted_retention_times[i] < 0.0) predicted_retention_times[i] = 0.0;
      else if (predicted_retention_times[i] > 1.0) predicted_retention_times[i] = 1.0;

      // predicted_retention_times[i] ->  needs scaling onto the column     
      SimCoordinateType rt_error = gsl_ran_gaussian(rnd_gen_, rt_shift_stddev) + rt_shift_mean;

      // shifting the retention time ensures that we do not have an elution profile that is
      // cut on the minima of the map
      SimCoordinateType retention_time = RTSimulation::gradient_front_offset_ + (predicted_retention_times[i] * (gradient_time_ - RTSimulation::gradient_total_offset_));

      std::cout << predicted_retention_times[i] << ", " << features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString() << ", abundance: " << features[i].getIntensity() << std::endl;
      
      features[i].setRT(retention_time + rt_error);
    } 
  }
  
  void RTSimulation::predictContaminantsRT(FeatureMapSim & contaminants)
  {
    // iterate of feature map
    for (Size i = 0; i < contaminants.size(); ++i)
    {
      // assign random retention time
      SimCoordinateType retention_time = gsl_ran_flat(rnd_gen_, 0, gradient_time_);
      contaminants[i].setRT(retention_time);
    }
  }
  
  void RTSimulation::updateMembers_()
  {
		rt_model_file_ = param_.getValue("rt_model_file");
		if (! File::readable( rt_model_file_ ) )
    { // look in OPENMS_DATA_PATH
      rt_model_file_ = File::find( rt_model_file_ );
    }
		gradient_time_ = param_.getValue("total_gradient_time");
    rt_sampling_rate_ = param_.getValue("sampling_rate");

    String column_preset = param_.getValue("column_condition:preset");
    if (column_preset == "poor")
    {
      distortion_    = 2.0;
      symmetry_down_ = -100;
      symmetry_up_   = +100;
    }
    else if (column_preset == "medium")
    {
      distortion_    = 1.0;
      symmetry_down_ = -60;
      symmetry_up_   = +60;
    }
    else if (column_preset == "good")
    {
      distortion_    = 0.0;
      symmetry_down_ = -15;
      symmetry_up_   = +15;
    }
    else 	// default is "none" so get user set parameters
    {
      distortion_    = param_.getValue("column_condition:distortion");
      symmetry_up_   = param_.getValue("column_condition:symmetry_up");
      symmetry_down_ = param_.getValue("column_condition:symmetry_down");
    }
    		
  }
  
  void RTSimulation::setDefaultParams_() 
  {
		defaults_.setValue("rt_column_on", "true", "Modelling of a rt column");
    defaults_.setValidStrings("rt_column_on", StringList::create("true,false"));
    
		// Column settings
    defaults_.setValue("total_gradient_time",2800.0,"the duration (in seconds) of the gradient");
    defaults_.setMinFloat("total_gradient_time", 800);
    // rt parameters
    defaults_.setValue("sampling_rate", 2.0, "Time interval (in seconds) between consecutive scans");

    defaults_.setValue("rt_model_file","examples/simulation/RTPredict.model","SVM model for retention time prediction");
    
    // rt error
    defaults_.setValue("rt_shift_mean",0,"Mean shift in retention time [s]");
    defaults_.setValue("rt_shift_stddev",50,"Standard deviation of shift in retention time [s]");     

    // column conditions
    defaults_.setValue("column_condition:preset","medium","LC condition (none|good|medium|poor) if set to none the explicit values will be used.");
    StringList valid_presets = StringList::create("none,good,medium,poor");
    defaults_.setValidStrings("column_condition:preset", valid_presets);

    defaults_.setValue("column_condition:distortion", 1.0, "LC distortion (used only if preset is set to 'none')");
    defaults_.setValue("column_condition:symmetry_up", -60.0, "LC symmetry up (used only if preset is set to 'none')");
    defaults_.setValue("column_condition:symmetry_down", +60.0, "LC symmetry down (used only if preset is set to 'none')");
    
		defaultsToParam_();
  }
  
  bool RTSimulation::isRTColumnOn() const
  {
    return (param_.getValue("rt_column_on") == "true");
  }
  
  SimCoordinateType RTSimulation::getGradientTime() const
  {
    return gradient_time_;
  }
  
  void RTSimulation::createExperiment_(MSSimExperiment & experiment)
  {
    std::cout << "create experiment ... ";
    experiment = MSSimExperiment();
    
    if (isRTColumnOn())
    {
      Size number_of_scans = Size(gradient_time_ / rt_sampling_rate_);
      experiment.resize(number_of_scans);

      DoubleReal current_scan_rt = rt_sampling_rate_;
      for(MSSimExperiment::iterator exp_it = experiment.begin();
          exp_it != experiment.end();
          ++exp_it)
      {
        (*exp_it).setRT(current_scan_rt);
        // dice & store distortion
        DoubleReal distortion = exp(gsl_ran_flat (rnd_gen_, -distortion_, +distortion_));
        (*exp_it).setMetaValue("distortion", distortion);
        
        // TODO: smooth?!
        
        /** WARNING: OLD CODE       
				// moving average filter (width 3), implemented very inefficiently (guess why!)
				if ( distortion_ != 0.0 ) // otherwise we want perfect EMG shape!
				{
					ElutionModel::ContainerType tmp;
					tmp.resize(data.size());
					for ( Size i = 1; i < data.size() - 1; ++i )
					{
						tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
					}
					for ( Size i = 1; i < data.size() - 1; ++i )
					{
						data[i] = tmp[i];
					}

					const int num_rounds = 10;
					for ( int rounds = 0; rounds < num_rounds; ++rounds)
					{
					//TODO: WTF!
						data.swap(tmp);
						for ( Size i = 1; i < data.size() - 1; ++i )
						{
							tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
						}
						for ( Size i = 1; i < data.size() - 1; ++i )
						{
							data[i] = tmp[i];
						}
					}

				}
				**/
                
        current_scan_rt += rt_sampling_rate_;
      }
    }
    else
    {
      experiment.resize(1);
      experiment[0].setRT(-1);
    }
    std::cout << "done\n";
  }
    
}
