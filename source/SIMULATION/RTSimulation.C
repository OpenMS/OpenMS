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
  
  void RTSimulation::no_rt_column(FeatureMapSim & features)
  {
    for(FeatureMapSim::iterator fIt = features.begin(); fIt != features.end();
        ++fIt)
    {
      (*fIt).setRT(-1);
    }
  }
  /**
   @brief Gets a feature map containing the peptides and predicts for those the retention times
   */
  void RTSimulation::predict_rt(FeatureMapSim & features)
  {
    Int doPredict = param_.getValue("rt_column_on");
    if(doPredict == 1)
    {
      svm_predict(features);
    }
    else
    {
      no_rt_column(features);
    }    
  }
  
  /**
   @brief Gets a feature map containing the peptides and predicts for those the retention times
   */
  void RTSimulation::svm_predict(FeatureMapSim & features)
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
    
    if (! File::readable( rtModelFile_ ) )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "RTSimulation got invalid parameter. 'rt_model_file' " + rtModelFile_ + " is not readable");
    }
    else
    {
      cout << "Predicting RT..    " << endl;
    }   
    
    // not that elegant...
    vector< String > peptidesVector(features.size());

    for (Size i = 0; i < features.size(); ++i)
    {
      peptidesVector[i] = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
    }
    
    svm.loadModel(rtModelFile_);
    
    // load additional parameters
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      String add_paramfile = rtModelFile_ + "_additional_parameters";
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
    String sample_file = rtModelFile_ + "_samples";
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
    
    for (UInt i = 0; i < peptidesVector.size(); i++)
    {
      // TODO: remove those peptides + debug out removed ones??
      if (predicted_retention_times[i] < 0.0) predicted_retention_times[i] = 0.0;
      else if (predicted_retention_times[i] > 1.0) predicted_retention_times[i] = 1.0;

      // predicted_retention_times[i] ->  needs scaling onto the column     
      SimCoordinateType rt_error = gsl_ran_gaussian(rnd_gen_, rt_shift_stddev) + rt_shift_mean;

      // shifting the retention time ensures that we do not have an elution profile that is
      // cut on the minima of the map
      SimCoordinateType retention_time = RTSimulation::gradient_front_offset_ + (predicted_retention_times[i] * (gradientTime_ - RTSimulation::gradient_total_offset_));

      std::cout << predicted_retention_times[i] << ", " << features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString() << ", abundance: " << features[i].getIntensity() << std::endl;
      
      features[i].setRT(retention_time + rt_error);
    } 
  }
  
  void RTSimulation::predict_contaminants_rt(FeatureMapSim & contaminants)
  {
    // iterate of feature map
    for (Size i = 0; i < contaminants.size(); ++i)
    {
      // assign random retention time
      SimCoordinateType retention_time = gsl_ran_flat(rnd_gen_, 0, gradientTime_);
      contaminants[i].setRT(retention_time);
    }
  }
  
  void RTSimulation::updateMembers_()
  {
    rtModelFile_ = param_.getValue("rt_model_file");
    gradientTime_ = param_.getValue("total_gradient_time");
  }
  
  void RTSimulation::setDefaultParams_() 
  {
    // we need to further integrate this .. currently it is ignored
    // TODO: we need to propagate this also to the product model which generates the signal
    // but how
    defaults_.setValue("rt_column_on",1,"Modelling of a rt column (0 = disabled, 1 = enabled)");
    
    // Column settings
    defaults_.setValue("total_gradient_time",2800.0,"the duration (in seconds) of the gradient");
    defaults_.setMinFloat("total_gradient_time", 800);
    defaults_.setValue("rt_model_file","file","SVM model for retention time prediction");
    
    // rt error
    defaults_.setValue("rt_shift_mean",0,"Mean shift in retention time [s]");
    defaults_.setValue("rt_shift_stddev",50,"Standard deviation of shift in retention time [s]");     

		defaultsToParam_();
  }
  
  bool RTSimulation::isRTColumnOn() const
  {
    Int isRTColumnOn = param_.getValue("rt_column_on");
    return (isRTColumnOn == 1);
  }
  
  SimCoordinateType RTSimulation::getGradientTime() const
  {
    return gradientTime_;
  }
}
