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
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"
#include "SqtSanityCheck.h"

// Class doing the sanity check for non SQT file condition
// In most sence a place holder as very little logic is build-in in this case

SanityCheck::SanityCheck() :
  initPositives_(0), pTestset_(NULL), pTrainset_(NULL), concatenatedSearch_(true),
  test_fdr_(0.01), initial_train_fdr_(0.01) {
}

SanityCheck::~SanityCheck() {
}

bool SanityCheck::overRule = false;
string SanityCheck::initWeightFN = "";
int SanityCheck::initDefaultDir = 0;
string SanityCheck::initDefaultDirName = "";
vector<double> SanityCheck::default_weights = vector<double>();

/**
 * Returns an instance of the appropriate sanity check based on information
 * contained in otherCall.
 */
SanityCheck* SanityCheck::initialize(string otherCall){
  if(initWeightFN != "" || initDefaultDirName.size() > 0) {
    return new SanityCheck();
  } else if (otherCall.find(SqtSanityCheck::fingerPrint)!= string::npos){
    return new SqtSanityCheck();
  } else {
    return new SanityCheck();
  }
}

void SanityCheck::checkAndSetDefaultDir() {
  if (!initDefaultDir && initDefaultDirName.size() > 0) {
    int sign = 1;
    if (initDefaultDirName[0] == '-') {
      initDefaultDirName.erase(0,1);
      sign = -1;
    }
    initDefaultDir = sign * DataSet::getFeatureNames().getFeatureNumber(initDefaultDirName);
    if (initDefaultDir == 0) {
      ostringstream temp;
      temp << "ERROR: Initial direction feature name \"" << initDefaultDirName << "\" not found" << std::endl;
      throw MyException(temp.str());
    }
  }
}

int SanityCheck::getInitDirection(vector<Scores>& testset,
                                  vector<Scores>& trainset,
                                  const Normalizer* pNorm,
                                  vector<vector<double> > & w,
                                  double test_fdr,
                                  double initial_train_fdr) {
  pTestset_ = &testset;
  pTrainset_ = &trainset;
  test_fdr_= test_fdr;
  initial_train_fdr_ = initial_train_fdr;
  if (initWeightFN.size() > 0) {
    vector<double> ww(FeatureNames::getNumFeatures() + 1);
    ifstream weightStream(initWeightFN.data(), ios::in);
    if (weightStream.is_open()) {
      readWeights(weightStream, ww);
      weightStream.close();
      assert(pNorm);
      pNorm->normalizeweight(ww, w[0]);
      for (size_t set = 1; set < w.size(); ++set) {
        copy(w[0].begin(), w[0].end(), w[set].begin());
      }
    } else {
      std::cerr << "WARNING: Could not find weights input file " << initWeightFN
                << ". Using default weights instead." << std::endl;
      getDefaultDirection(w);
    }
  } else {
    getDefaultDirection(w);
  }
  initPositives_ = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    initPositives_ += (*pTestset_)[set].calcScoresAndQvals(w[set], test_fdr);
  }
  return initPositives_;
}

void SanityCheck::getDefaultDirection(vector<vector<double> >& w) {
    
  //If I have not been given a initial direction
  if (!initDefaultDir) {
    //if the default_weights from pin.xml are not present
    if (default_weights.size() == 0) {
      // Set init direction to be the most discriminative direction
      for (size_t set = 0; set < w.size(); ++set) {
        if (VERB > 1) std::cerr << "Split " << set + 1 << ":\t";
        calcInitDirection(w[set], set);
      }
    } else {
      // I want to assign the default vector that is present in the input file
      for (size_t set = 0; set < w.size(); ++set) {
        for (size_t ix = 0; ix < w[set].size(); ix++) {
          w[set][ix] = 0;
          if (ix < default_weights.size()){
              w[set][ix] = default_weights[ix];
          }          
        }
      }
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      for (size_t ix = 0; ix < w[set].size(); ix++) {
        w[set][ix] = 0;
      }
      w[set][static_cast<std::size_t>(abs(initDefaultDir) - 1)] = (initDefaultDir < 0 ? -1 : 1);
    }
  }
}

void SanityCheck::calcInitDirection(vector<double>& wSet, size_t set) {
  (*pTrainset_)[set].getInitDirection(initial_train_fdr_, wSet);
}

bool SanityCheck::validateDirection(vector<vector<double> >& w) {
  if (!pTestset_) {
    cerr << "SanityCheck wrongly configured" << endl;
    return false;
  }
  int overFDR = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    overFDR += (*pTestset_)[set].calcScoresAndQvals(w[set], test_fdr_);
  }
  if (overFDR <= 0) {
    cerr << "No targets found with q<" << test_fdr_ << endl;
    resetDirection(w);
    return false;
  }
  if (initPositives_ > overFDR) {
    cerr << "Less identifications (" << overFDR << " vs " << initPositives_
        << ") after percolator processing than before processing" << endl;
    resetDirection(w);
    return false;
  }
  if (initDefaultDir) {
    for (size_t set = 0; set < w.size(); ++set) {
      if (w[set][static_cast<std::size_t>(abs(initDefaultDir) - 1)] * initDefaultDir <= 0) {
        if (VERB > 1) {
          cerr
            << "Warning, wrong sign of the weight for main scoring direction"
            << endl;
        }
        resetDirection(w);
        return false;
      }
    }
  }
  return true;
}

void SanityCheck::readWeights(istream& weightStream, vector<double>& w) {
  std::string line;
  // Find first line with weights, but skip this and get second line containing raw features
  while (std::getline(weightStream, line)) {
    if (line[0] == '-' || (line[0] >= '0' && line[0] <= '9')) break;
  }
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream >> w[ix];
  }
  cerr << "Read weights from file" << endl;
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    cerr << w[ix] << "\t";
  }
  cerr << endl;
}

void SanityCheck::resetDirection(vector<vector<double> >& w) {
  if (!overRule) {
    cerr << "Resetting score vector, using default vector. Use --override flag to prevent this." << endl;
    getDefaultDirection(w);
  }
}

