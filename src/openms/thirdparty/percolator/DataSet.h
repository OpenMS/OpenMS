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
#include <string_view>
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
#include "Globals.h"
#include "PSMDescription.h"
#include "PsmStore.h"
#include "TabReader.h"
#include "FeatureNames.h"
#include "FeatureMemoryPool.h"

namespace OpenMS { namespace Internal { namespace Percolator {



// Optional columns in tab delimited input
enum OptionalField {
  SCANNR, EXPMASS, CALCMASS, RETTIME, DELTAMASS, FILENAME
};

// The fields of a PSM line before its features, see DataSet::readPsmHeader
struct PsmHeader {
  // the scan number, the experimental mass and the line number
  PSMDescription psm;
  std::string_view id;  // a view into the line
  int label = 0;
  // the starts of the fields that DataSet::readPsmBody parses, if present
  const char* calcMassField = NULL;
  const char* retentionTimeField = NULL;
  bool hasFileName = false;
  std::string_view fileName;
};

class DataSet {
 public:
  DataSet();
  virtual ~DataSet();
  DataSet(const DataSet&) = delete;
  DataSet& operator=(const DataSet&) = delete;
  
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
  

  void fillScores(std::vector<ScoreHolder>& scores);
  void fillFeatures(std::vector<double*>& features);
  void clearFeaturePointers();
  
  // Reads a PSM line into a new PSM of the store, without checking the
  // fields before the features, which SetHandler::getScanId does. For the
  // lines that the fast path below refuses.
  static LabelType readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool,
    PsmStore& psmStore, const std::string& decoyPrefix);

  // A fast path for readPsm that reads the fields before the features only
  // once; see DataSet.cpp. Both return false if the line has to be read with
  // SetHandler::getScanId and readPsm instead.
  static bool readPsmHeader(TabReader& reader, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, PsmHeader& header);
  static bool readPsmBody(TabReader& reader, PsmHeader& header,
    bool readProteins, PSMDescription*& myPsm, LabelType& label,
    FeatureMemoryPool& featurePool, PsmStore& psmStore,
    const std::string& decoyPrefix);
  
  void registerPsm(PSMDescription* myPsm);

  // lets the warnings that are printed once per run be printed again
  static void resetWarnings() {
    decoyWarningTripped_ = false;
    flankWarningTripped_ = false;
  }
  
 protected:   
  
  std::vector<PSMDescription*> psms_;  // in input order, owned by a PsmStore
  LabelType label_ = LabelType::UNDEFINED;
  
  static void checkUnflankedPeptide(std::string_view peptide,
                                    std::string_view psmId,
                                    unsigned int lineNr);

  static FeatureNames featureNames_;
  static bool decoyWarningTripped_;
  static bool flankWarningTripped_;

  static void readFeaturesAndProteins(TabReader& reader, LabelType label,
    bool readProteins, PSMDescription& psm, std::string_view id,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool,
    PsmStore& psmStore, const std::string& decoyPrefix);
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*DATASET_H_*/
