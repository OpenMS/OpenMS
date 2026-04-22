/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef DATASET_H_
#define DATASET_H_

#include <string>
#include <cassert>
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <cerrno>
#include <random>

#include "Scores.h"
#include "ResultHolder.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "FeatureNames.h"
#include "FeatureMemoryPool.h"
#include "ProteinProbEstimator.h"
#include "TabFileValidator.h"


// Optional columns in tab delimited input
enum OptionalField {
  SCANNR, EXPMASS, CALCMASS, RETTIME, DELTAMASS, FILENAME
};

class DataSet {
 public:
  DataSet();
  virtual ~DataSet();
  
  void initFeatureTables(const unsigned int numFeatures);
  
  void inline setLabel(LabelType l) { label_ = l; }
  LabelType inline getLabel() const { return label_; }
  
  unsigned int inline getSize() const { return static_cast<unsigned int>(psms_.size()); }
    
  static FeatureNames& getFeatureNames() { return featureNames_; }
  static void resetFeatureNames() { 
    featureNames_ = FeatureNames();
    FeatureNames::resetNumFeatures();
  }
  static unsigned getNumFeatures() { return static_cast<unsigned>(featureNames_.getNumFeatures()); }
  
  bool writeTabData(std::ofstream& out);
  
  void print_10features();
  void print_features();

  void fillScores(std::vector<ScoreHolder>& scores);
  void fillFeatures(std::vector<double*>& features);
  
  void readPsm(const std::string& line, const unsigned int lineNr,
               const std::vector<OptionalField>& optionalFields, 
               FeatureMemoryPool& featurePool, std::string decoyPrefix);
  static LabelType readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool, std::string decoyPrefix);
  
  void registerPsm(PSMDescription* myPsm);
  
 protected:   
  
  std::vector<PSMDescription*> psms_;
  LabelType label_;
  
  static FeatureNames featureNames_;
  static bool decoyWarningTripped_;
};

#endif /*DATASET_H_*/
