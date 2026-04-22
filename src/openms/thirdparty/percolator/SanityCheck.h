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
#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_

class Scores;
class Normalizer;

class SanityCheck {
 public:
  SanityCheck();
  virtual ~SanityCheck();
  static SanityCheck* initialize(string otherCall);
  void readWeights(istream& weightStream, vector<double>& w);
  void checkAndSetDefaultDir();
  int getInitDirection(vector<Scores>& testset,
                       vector<Scores> &trainset, const Normalizer* pNorm,
                       vector<vector<double> >& w, double test_fdr,
                       double initial_train_fdr);
  
  virtual bool validateDirection(vector<vector<double> >& w);
  void resetDirection(vector<vector<double> >& w);

  static void setInitWeightFN(string fn) {
    initWeightFN = fn;
  }
  static void setInitDefaultDir(int dir) {
    initDefaultDir = dir;
  }
  static void setInitDefaultDirName(const string & dirName) {
    initDefaultDirName = dirName;
  }
  static void setOverrule(bool orl) {
    overRule = orl;
  }
  
  static void addDefaultWeights(vector<double> __default_weights) {
      default_weights = __default_weights;
  }
  
  static vector<double>& getDefaultWeights() {
      return default_weights;
  }
  
  bool concatenatedSearch() const {
    return concatenatedSearch_;
  }
  
  void setConcatenatedSearch(const bool s) {
    concatenatedSearch_ = s;
  }
  
 protected:
  virtual void getDefaultDirection(vector<vector<double> >& w);
  virtual void calcInitDirection(vector<double>& wSet, size_t set);
  
  int initPositives_;
  double test_fdr_, initial_train_fdr_;
  
  static bool overRule;
  static string initWeightFN;
  static string initDefaultDirName;
  static int initDefaultDir; // Default Direction, 0=do not use,
  // positive integer = feature number,
  // negative integer = lower score better
  
  vector<Scores> *pTestset_, *pTrainset_;
  
  //this vector is meant to have the initial directions givn in the pin.xml
  static vector<double> default_weights;
  
  // input from concatenated search, i.e. not separate target and decoy searches
  bool concatenatedSearch_;
};

#endif /*SANITYCHECK_H_*/
