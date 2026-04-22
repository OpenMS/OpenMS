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
#ifndef SCORE_HOLDER_H_
#define SCORE_HOLDER_H_

#include "LabelType.h"
#include "PSMDescription.h"

class Scores;  // forward declaration

/*
 * ScoreHolder is a class that provides a way to assign score value to a
 * PSMDescription and have a way to compare PSMs based on the assigned
 * score value and output them into the stream.
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 *
 */
class ScoreHolder {
 public:
  double score, q, pep, p;
  PSMDescription* pPSM;
  LabelType label;

  ScoreHolder()
      : score(0.0),
        q(0.0),
        pep(0.0),
        p(0.0),
        label(LabelType::UNDEFINED),
        pPSM(NULL) {}
  ScoreHolder(const double s, const LabelType l, PSMDescription* psm = NULL)
      : score(s), q(0.0), pep(0.0), p(0.0), label(l), pPSM(psm) {}
  virtual ~ScoreHolder() {}

  std::pair<double, bool> toPair() const {
    return pair<double, bool>(score, isTarget());
  }

  inline bool isTarget() const {
    return label == LabelType::TARGET || label == LabelType::PSEUDO_TARGET;
  }
  inline bool isDecoy() const { return label == LabelType::DECOY; }
  inline PSMDescription* getPSM() const { return pPSM; }
  void printPSM(ostream& os, bool printDecoys, bool printExpMass);
  void printPepXML(ostream& os, map<char, float>& aaWeight, int index);
  void printPeptide(ostream& os,
                    bool printDecoys,
                    bool printExpMass,
                    Scores& fullset);
  std::string getCharge(std::string id);
};
struct lessThanBaseName {
  inline bool operator()(const ScoreHolder& struct1,
                         const ScoreHolder& struct2) {
    std::string id1 = struct1.pPSM->getId();
    std::string id2 = struct2.pPSM->getId();
    return (id1 < id2);
  }
};
bool operator>(const ScoreHolder& one, const ScoreHolder& other);
bool operator<(const ScoreHolder& one, const ScoreHolder& other);

struct lexicOrderProb {
  static int compStrIt(const std::string& sv1, const std::string& sv2) {
    auto it1 = sv1.begin(), end1 = sv1.end();
    auto it2 = sv2.begin(), end2 = sv2.end();

    for (; it1 != end1 && it2 != end2; ++it1, ++it2) {
      if (*it1 < *it2)
        return 1;
      if (*it2 < *it1)
        return -1;
    }
    if (it2 != end2)
      return 1;
    if (it1 != end1)
      return -1;
    return 0;
  }

  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    const std::string& fullSeqX = x.pPSM->getFullPeptideSequence();
    const std::string& fullSeqY = y.pPSM->getFullPeptideSequence();

    std::string peptSeqX = fullSeqX.substr(2, fullSeqX.size() - 4);
    std::string peptSeqY = fullSeqY.substr(2, fullSeqY.size() - 4);

    int peptCmp = compStrIt(peptSeqX, peptSeqY);
    return ((peptCmp == 1) || ((peptCmp == 0) && (x.label > y.label)) ||
            ((peptCmp == 0) && (x.label == y.label) && (x.score > y.score)));
  }
};

/**
 * Computes a fast hash for unsigned int using the function suggested at:
 * https://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
 * Note that std::hash is sometimes implemented as the identity function,
 * defeating the purpose of random shuffling.
 */
inline unsigned int fast_uint_hash(unsigned int x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x;
}

/**
 * Orders ScoreHolders by a hash computed from the specFileNr and scan.
 * Uses the hash combination function h1 ^ (h2 << 1) suggested in
 * https://en.cppreference.com/w/cpp/utility/hash
 */
struct OrderScanHash {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    size_t hash_x = fast_uint_hash(x.pPSM->specFileNr) ^
                    (fast_uint_hash(x.pPSM->scan) << 1);
    size_t hash_y = fast_uint_hash(y.pPSM->specFileNr) ^
                    (fast_uint_hash(y.pPSM->scan) << 1);
    return hash_x < hash_y;
  }
};

struct OrderScanMassCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    if (x.pPSM->specFileNr != y.pPSM->specFileNr)
      return x.pPSM->specFileNr < y.pPSM->specFileNr;
    if (x.pPSM->scan != y.pPSM->scan)
      return x.pPSM->scan < y.pPSM->scan;
    if (x.pPSM->expMass != y.pPSM->expMass)
      return x.pPSM->expMass < y.pPSM->expMass;
    return x.score > y.score;
  }
};

struct OrderScanMassLabelCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    if (x.pPSM->specFileNr != y.pPSM->specFileNr)
      return x.pPSM->specFileNr < y.pPSM->specFileNr;
    if (x.pPSM->scan != y.pPSM->scan)
      return x.pPSM->scan < y.pPSM->scan;
    if (x.pPSM->expMass != y.pPSM->expMass)
      return x.pPSM->expMass < y.pPSM->expMass;
    if (x.label != y.label)
      return x.label > y.label;
    return x.score > y.score;
  }
};

struct OrderScanLabel {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    if (x.pPSM->specFileNr != y.pPSM->specFileNr)
      return x.pPSM->specFileNr < y.pPSM->specFileNr;
    if (x.pPSM->scan != y.pPSM->scan)
      return x.pPSM->scan < y.pPSM->scan;
    return x.label > y.label;
  }
};

struct UniqueScanMassCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    return (x.pPSM->specFileNr == y.pPSM->specFileNr) &&
           (x.pPSM->scan == y.pPSM->scan) &&
           (x.pPSM->expMass == y.pPSM->expMass);
  }
};

struct UniqueScanMassLabelCharge {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    return (x.pPSM->specFileNr == y.pPSM->specFileNr) &&
           (x.pPSM->scan == y.pPSM->scan) && (x.label == y.label) &&
           (x.pPSM->expMass == y.pPSM->expMass);
  }
};

struct UniqueScanLabel {
  bool operator()(const ScoreHolder& x, const ScoreHolder& y) const {
    return (x.pPSM->specFileNr == y.pPSM->specFileNr) &&
           (x.pPSM->scan == y.pPSM->scan) && (x.label == y.label);
  }
};

#endif /*SCORE_HOLDER_H_*/