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

#include<OpenMS/SIMULATION/DetectabilitySimulation.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>


#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace OpenMS {

  DetectabilitySimulation::DetectabilitySimulation()
    : DefaultParamHandler("DetectabilitySimulation")
  {
    setDefaultParams_();
  }

  DetectabilitySimulation::DetectabilitySimulation(const DetectabilitySimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    updateMembers_(); 
  }

  DetectabilitySimulation& DetectabilitySimulation::operator = (const DetectabilitySimulation& source)
  {
    setParameters( source.getParameters() );
    updateMembers_();
    return *this;
  }
  
  DetectabilitySimulation::~DetectabilitySimulation()
  {}
  
  void DetectabilitySimulation::filterDetectability(FeatureMapSim & features)
  {
    if (param_.getValue("dt_simulation_on") == "true")
    {
      svm_filter(features);
    }
    else
    {
      no_filter(features);
    }
  }
  
  void DetectabilitySimulation::no_filter(FeatureMapSim & features)
  {  
    // set detectibility to 1.0 for all given peptides
    DoubleReal defaultDetectibility = 1.0;
    
    for(FeatureMapSim::iterator feature_it = features.begin();
        feature_it != features.end();
        ++feature_it) 
    {
      (*feature_it).setMetaValue("detectibility", defaultDetectibility );
    }     
  }
    
  void DetectabilitySimulation::svm_filter(FeatureMapSim & features)
  {
    
    // The support vector machine
		SVMWrapper svm_;

    // initialize support vector machine
    LibSVMEncoder encoder;
    svm_problem* training_data = NULL;
    UInt k_mer_length = 0;
    DoubleReal sigma = 0.0;
    UInt border_length = 0;
    
		if (File::readable(dtModelFile_))
		{
    	svm_.loadModel(dtModelFile_);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "DetectibilitySimulation got invalid parameter. 'dt_model_file' " + dtModelFile_ + " is not readable");
    }
		
		// load additional parameters
		if (svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      String add_paramfile = dtModelFile_ + "_additional_parameters";
      if (! File::readable( add_paramfile ) )
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "DetectibilitySimulation: SVM parameter file " + add_paramfile + " is not readable");
      }
      
      Param additional_parameters;
      additional_parameters.load(add_paramfile);
      
      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
          && svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "DetectibilitySimulation: No border length defined in additional parameters file.");
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
          && svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "DetectibilitySimulation: No k-mer length defined in additional parameters file.");
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
      
      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
          && svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "DetectibilitySimulation: No sigma defined in additional parameters file.");
      }
      
      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
    }
    
		if (File::readable(dtModelFile_))
		{
    	svm_.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
    	svm_.setParameter(SVMWrapper::SIGMA, sigma);
    	// to obtain probabilities
    	svm_.setParameter(SVMWrapper::PROBABILITY, 1);
		}
    // loading training data
    String sample_file = dtModelFile_ + "_samples";
    if (File::readable(sample_file))
    {
    	training_data = encoder.loadLibSVMProblem(sample_file);
    	svm_.setTrainingSample(training_data);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "DetectibilitySimulation: SVM sample file " + sample_file + " is not readable");
    }
    // transform featuremap to peptides vector
    vector< String > peptides_vector(features.size());
    for(Size i = 0; i < features.size(); ++i)
    {
      peptides_vector[i] = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
    }
    
    
    cout << "Predicting peptide detectabilities..    " << endl;
    
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    
    // Encoding test data
    vector<DoubleReal> probs;
    probs.resize(peptides_vector.size(), 0);

    svm_problem* prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptides_vector, probs,
                                                                                     k_mer_length,
                                                                                     allowed_amino_acid_characters,
                                                                                     svm_.getIntParameter(SVMWrapper::BORDER_LENGTH));
    
    vector<DoubleReal> labels;
    vector<DoubleReal> detectabilities;
    svm_.getSVCProbabilities(prediction_data, detectabilities, labels);
    
    cout << "Done." << endl;
    
    delete prediction_data;
    
#ifdef DEBUG_SIM
    cout << "----------------------------------------------------------------" << endl;
    cout << "Predicted detectabilities:" << endl;
#endif
    
    // copy all meta data stored in the feature map
    FeatureMapSim temp_copy(features); 
    temp_copy.clear();
    
    for (Size i = 0; i < peptides_vector.size(); ++i)
    {

      if (detectabilities[i] > min_detect_)
      {
        features[i].setMetaValue("detectability", detectabilities[i] );
        temp_copy.push_back(features[i]);
      }
#ifdef DEBUG_SIM
      cout << detectabilities[i] << " " << min_detect_ << endl;
#endif
    } 
    
    features.swap(temp_copy);
  }

  void DetectabilitySimulation::setDefaultParams_() 
  {
		defaults_.setValue("dt_simulation_on", "false", "Modelling detectibility");
    defaults_.setValidStrings("dt_simulation_on", StringList::create("true,false"));
    defaults_.setValue("min_detect",0.5,"Minimum peptide detectability accepted");
    defaults_.setValue("dt_model_file","","SVM model for peptide detectability prediction");
    defaultsToParam_();
  }
  
  void DetectabilitySimulation::updateMembers_()
  {
    min_detect_ = param_.getValue("min_detect");
    dtModelFile_ = param_.getValue("dt_model_file");
  }  
}
