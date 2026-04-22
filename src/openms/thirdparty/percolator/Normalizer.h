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
#ifndef NORMALIZER_H_
#define NORMALIZER_H_
#include <set>
#include <vector>
#include <iostream>

using namespace std;

class Normalizer {
 public:
  virtual ~Normalizer();
  virtual void setSet(vector<double*>& /* featuresV */,
                      vector<double*>& /* rtFeaturesV */, size_t /* numFeatures */,
                      size_t /* numRetentionFeatures*/ ) {}
  virtual void updateSet(vector<double*>& /* featuresV */, size_t /* offset */,
                         size_t /* numFeatures */) {}
  
  void normalizeSet(vector<double*>& featuresV,
                    vector<double*>& rtFeaturesV);
  void normalizeSet(vector<double*>& featuresV,
                    size_t offset, size_t numFeatures);
  void normalize(const double* in, double* out, size_t offset,
                 size_t numFeatures);
  inline double normalize(const double in, size_t index) {
    return (in - sub[index]) / div[index];
  }
  virtual void unnormalizeweight(const vector<double>& /* in */,
                                 vector<double>& /* out */) const {}
  virtual void normalizeweight(const vector<double>& /* in */,
                               vector<double>& /* out */) const {}
  static void endScoreNormalizeWeights(const std::vector<double>& in, 
    std::vector<double>& out, double subScore, double scale);

  static Normalizer* getNormalizer();
  static void resetNormalizer() {
    theNormalizer = NULL;
  }
  static void setType(int type);
  const static int UNI = 0;
  const static int STDV = 1;
  const static int NONORM = 2;
  
  void resizeVecs(size_t size) {
    sub.resize(size, 0.0);
    div.resize(size, 1.0);
  }
  void setNumberRetentionFeatures(size_t numRF) {
    numRetentionFeatures = numRF;
  }
  void setNumFeatures(const size_t nf) {
    numFeatures = nf;
  }
  double* getSub() {
    return &sub[0];
  }
  double* getDiv() {
    return &div[0];
  }
  size_t* getNumRetFeatures() {
    return &numRetentionFeatures;
  }
  void printNumRetFeatures() {
    cout << "There are " << numRetentionFeatures << endl;
  }
  void printSub() {
    for (unsigned int i = 0; i < numRetentionFeatures; i++) {
      cout << sub[i] << " ";
    }
    cout << endl;
  }
  void printDiv() {
    for (unsigned int i = 0; i < numRetentionFeatures; i++) {
      cout << div[i] << " ";
    }
    cout << endl;
  }
  void SetSubDiv(const vector<double> s, const vector<double> d) {
    sub = s;
    div = d;
    numFeatures = 0;
    numRetentionFeatures = s.size();
  }
  vector<double> GetVSub() const { return sub; }
  vector<double> GetVDiv() const { return div; }
 protected:
  Normalizer();
  static Normalizer* theNormalizer;
  static int subclass_type;
  size_t numFeatures, numRetentionFeatures;
  vector<double> sub;
  vector<double> div;
};

#endif /*NORMALIZER_H_*/
