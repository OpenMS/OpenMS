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
// $Authors: Stephan Aiche Chris Bielow$
// --------------------------------------------------------------------------

#include<OpenMS/SIMULATION/DetectibilitySimulation.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>


#include <vector>
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

namespace OpenMS {

  DetectibilitySimulation::DetectibilitySimulation()
    : DefaultParamHandler("DetectibilitySimulation")
  {
    setDefaultParams_();
  }

  DetectibilitySimulation::DetectibilitySimulation(const DetectibilitySimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    updateMembers_(); 
  }

  DetectibilitySimulation& DetectibilitySimulation::operator = (const DetectibilitySimulation& source)
  {
    setParameters( source.getParameters() );
    updateMembers_();
    return *this;
  }
  
  DetectibilitySimulation::~DetectibilitySimulation()
  {}
  
  void DetectibilitySimulation::filterDetectibility(FeatureMap< > & features)
  {
    Int is_filter_active = param_.getValue("dt_simulation_on");
    if(is_filter_active == 1)
    {
      svm_filter(features);
    }
    else
    {
      no_filter(features);
    }
  }
  
  void DetectibilitySimulation::no_filter(FeatureMap< > & features)
  {  
    // set detectibility to 1.0 for all given peptides
    for (Size i = 0; i < features.size(); ++i)
    {
      features[i].setMetaValue("detectibility", 1.0 );
    }     
  }
    
  void DetectibilitySimulation::svm_filter(FeatureMap< > & features)
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
    	//cout << "Loading svm model: " << DtModelFile_ << endl;
    	svm_.loadModel(dtModelFile_);
			//cout << "Done. " << endl;
    }
		
		// load additional parameters
		if (svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      String add_paramfile = dtModelFile_ + "_additional_parameters";
      if (! File::readable( add_paramfile ) )
      {
        cout << "SVM parameter file " << dtModelFile_ << " not found or not readable" << endl;
        cout << "Not performing detectability prediction!" << endl;
      }
      
      Param additional_parameters;
      additional_parameters.load(add_paramfile);
      
      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
          && svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        cout << "No border length defined in additional parameters file." << endl;
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
          && svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        cout << "No k-mer length defined in additional parameters file. Aborting detectability prediction!" << endl;
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
      
      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
          && svm_.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        cout << "No sigma defined in additional parameters file. Aborting detectability prediction!" << endl;
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
    
    // transform featuremap to peptides vector
    vector< String > peptidesVector(features.size());
    for(Size i = 0; i < features.size(); ++i)
    {
      peptidesVector[i] = features[i].getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
    }
    
    
    cout << "Predicting peptide detectabilities..    " << endl;
    
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    
    // Encoding test data
    vector<DoubleReal> probs;
    probs.resize(peptidesVector.size(), 0);

    svm_problem* prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptidesVector, probs,
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
    
    FeatureMap< > tempCopy; 
    for (Size i = 0; i < peptidesVector.size(); ++i)
    {

      if (detectabilities[i] > min_detect_)
      {
        features[i].setMetaValue("detectibility", detectabilities[i] );
        tempCopy.push_back(features[i]);
      }
#ifdef DEBUG_SIM
      cout << detectabilities[i] << " " << min_detect_ << endl;
#endif
    } 
    
    features.swap(tempCopy);
  }

  void DetectibilitySimulation::setDefaultParams_() 
  {
    defaults_.setValue("min_detect",0.5,"minimum peptide detectability accepted");
    defaults_.setValue("dt_model_file","<file>","SVM model for peptide detectability prediction");
    defaults_.setValue("dt_simulation_on",1, "Modelling detectibility (0 = disabled, 1 = enabled)");
    defaultsToParam_();
  }
  
  void DetectibilitySimulation::updateMembers_()
  {
    min_detect_ = param_.getValue("min_detect");
    dtModelFile_ = param_.getValue("dt_model_file");
  }  
}
