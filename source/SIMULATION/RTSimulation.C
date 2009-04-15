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

#include<OpenMS/SIMULATION/RTSimulation.h>
#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <vector>
#include <iostream>
#include <utility>

using std::vector;
using std::cout;
using std::endl;
using std::make_pair;

namespace OpenMS {

  RTSimulation::RTSimulation()
    : DefaultParamHandler("RTSimulation")
  {
    // 
    defaults_.setValue("rt_column_on",1,"Modelling of a rt column (0 = disabled, 1 = enabled)");

    // Column settings
    defaults_.setValue("total_gradient_time",2000.0,"the duration (in seconds) of the gradient");
    defaults_.setValue("rt_sampling",2.0,"Time interval (in seconds) between consecutive scans");
    defaults_.setValue("rt_model_file","","");
    
		defaultsToParam_();    
  }

  RTSimulation::RTSimulation(const RTSimulation& source)
    : DefaultParamHandler(source)
  {}

  RTSimulation& RTSimulation::operator = (const RTSimulation& source)
  {
    return *this;
  }
  
  RTSimulation::~RTSimulation()
  {}
  
  /**
   @brief Gets a feature map containing the peptides and predicts for those the retention times
   */
  void RTSimulation::predict_rt(FeatureMap< > & input_features, FeatureMap< > & rt_features)
  {
    // TODO: maybe we should operate on the real data
    // TODO: compare results of this variant with the old predictions, just to ensure that everything works
    // TODO: define a MetaInfoField for sequence ..
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    SVMWrapper svm;
    LibSVMEncoder encoder;
    vector<DoubleReal> predicted_retention_times;
    
    svm_problem* training_data = NULL;
    svm_problem* prediction_data = NULL;
    
    UInt k_mer_length = 0;
    DoubleReal sigma = 0.0;
    UInt border_length = 0;
    
    cout << "Predicting RT..    " << endl;
    
    // not that elegant...
    vector< String > peptidesVector(input_features.size());

    //for(PeptideSequences::const_iterator seq_it = sample_.getPeptideSequences().begin();
    //    seq_it != sample_.getPeptideSequences().end();
    //    ++seq_it)
    for (Size i = 0; i < input_features.size(); ++i)
    {
      peptidesVector[i] = input_features[i].getMetaValue("sequence").toString();
      //peptidesVector.push_back(seq_it->first);
    }
    
    svm.loadModel(RTModelFile_);
    
    // load additional parameters
    if (svm.getIntParameter(KERNEL_TYPE) == OLIGO)
    {
      String add_paramfile = RTModelFile_ + "_additional_parameters";
      if (! File::readable( add_paramfile ) )
      {
        cout << "SVM parameter file " << add_paramfile << " not found or not readable" << endl;
        cout << "Aborting RT prediction!" << endl;
        //TODO where is the "return;" ?
        //TODO why is rtTable not filled with default values first?! (empty rtTable leads to empty map?! -->  see LCMSSim::run())
				// Ole: these are valid suggestions, but they should be adressed when the simulator is modified to 
				// incorporate MS/MS spectra.
      }
      
      Param additional_parameters;
      additional_parameters.load(add_paramfile);
      
      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
          && svm.getIntParameter(KERNEL_TYPE) == OLIGO)
      {
        cout << "No border length defined in additional parameters file. Aborting RT prediction!" << endl;
        return;
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
          && svm.getIntParameter(KERNEL_TYPE) == OLIGO)
      {
        cout << "No k-mer length defined in additional parameters file. Aborting RT prediction!" << endl;
        return;
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();
      
      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
          && svm.getIntParameter(KERNEL_TYPE) == OLIGO)
      {
        cout << "No sigma defined in additional parameters file. Aborting RT prediction!" << endl;
        return;
      }
      
      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
    }
    
    svm.setParameter(BORDER_LENGTH, (Int) border_length);
    svm.setParameter(SIGMA, sigma);
    
    // Encoding test data
    vector<DoubleReal> rts;
    rts.resize(peptidesVector.size(), 0);
    prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptidesVector, rts, k_mer_length, allowed_amino_acid_characters, border_length);
    
    // loading training data
    String sample_file = RTModelFile_ + "_samples";
    if (! File::readable( sample_file ) )
    {
      cout << "SVM sample file " << sample_file << " not found or not readable" << endl;
      cout << "Aborting RT prediction!" << endl;
    }
    training_data = encoder.loadLibSVMProblem(sample_file);
    
    svm.setTrainingSample(training_data);
    svm.predict(prediction_data, predicted_retention_times);
    
    cout << "Done." << endl;
    delete training_data;
    delete prediction_data;
    
    rt_features = input_features;
    for (UInt i = 0; i < peptidesVector.size(); i++)
    {
      if (predicted_retention_times[i] < 0.0) predicted_retention_times[i] = 0.0;
      else if (predicted_retention_times[i] > 1.0) predicted_retention_times[i] = 1.0;

      rt_features[i].setMetaValue("rt_time",predicted_retention_times[i]);
    }
    // Create the retention time table:
    // RT : pointer to peptide sequence : rel. abundance
    /*
    for (UInt i = 0; i < peptidesVector.size(); i++)
    {
      PeptideSequences::const_iterator pit = sample_.getPeptideSequences().find( peptidesVector[i] );
      
      // check for rt outlier and rectify
      if (predicted_retention_times[i] < 0.0) predicted_retention_times[i] = 0.0;
      else if (predicted_retention_times[i] > 1.0) predicted_retention_times[i] = 1.0;
      
      rtTable.insert( make_pair( predicted_retention_times[i], pit ) );
    }
    */
#ifdef DEBUG_SIM
    cout << "---------------------------------------------------------" << endl;
    cout << "Content of retention time table:" << endl;
    for (RTTable::const_iterator rtTableRow = rtTable.begin();
         rtTableRow != rtTable.end();
         ++rtTableRow)
    {
      cout << rtTableRow->first << ", ";
      cout << ( (rtTableRow->second)->first ) << " , abundance: " << ( (rtTableRow->second)->second ) << endl;
    }
    cout << "---------------------------------------------------------" << endl;
#endif  
  
  }
}
