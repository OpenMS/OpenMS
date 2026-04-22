
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

#include "FeatureNames.h"
#include "Globals.h"

size_t FeatureNames::numFeatures = 0;

FeatureNames::FeatureNames() {
  minCharge = 100;
  maxCharge = -1;
  chargeFeatNum = -1;
  enzFeatNum = -1;
  numSPFeatNum = -1;
  ptmFeatNum = -1;
  intraSetFeatNum = -1;
  quadraticFeatNum = -1;
}

FeatureNames::~FeatureNames() {
}

void FeatureNames::initFeatures() {
  setNumFeatures(featureNames.size());
  if (VERB>2) {
    std::cerr << "in FeatureNames::initFeatures\n";
  }
  if (VERB>1) {
    std::cerr << "Features:\n";
    std::copy( featureNames.begin(), featureNames.end(),
        std::ostream_iterator<std::string>(std::cerr, " "));
    std::cerr << "\n";
  }
  if (VERB>2) {
    std::cerr << "end of FeatureNames::initFeatures\n";
  }
}

string FeatureNames::getFeatureNames() {
  int n = (int)featureNames.size();
  ostringstream oss;
  if (!featureNames.empty()) {
    std::size_t featNum = 0;
    oss << featureNames[featNum++];
    for (; static_cast<int>(featNum) < n; ++featNum) {
      oss << "\t" << featureNames[featNum];
    }
  }
  return oss.str();
}

int FeatureNames::getFeatureNumber(const string& featureName) {
  for (unsigned int fnum = 0; fnum < featureNames.size(); ++fnum) {
    // there is no easy case insensitive string compare in c++, so go by char
    if (featureNames[fnum].size() == featureName.size()) {
      bool isEqual = true;
      for (unsigned int i = 0; i < featureName.size(); ++i) {
        if (std::tolower(featureName[i]) != std::tolower(featureNames[fnum][i])) {
          isEqual = false;
          break;
        }
      }
      if (isEqual) return static_cast<int>(fnum + 1);
    } 
  }
  return 0;
}
