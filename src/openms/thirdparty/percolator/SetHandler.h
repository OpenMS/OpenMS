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

#ifndef SETHANDLER_H_
#define SETHANDLER_H_

#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include <string>
#include <cctype>
#include <locale>
#include <queue>
#include <climits>

#include "DataSet.h"
#include "Normalizer.h"
#include "Scores.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "PseudoRandom.h"
#include "FeatureMemoryPool.h"
#include "ScanIdTable.h"

namespace OpenMS { namespace Internal { namespace Percolator {

// OPENMS_MODIFIED: the input files of a multi-file input, from ValidateTabFile.h (the multi-file reader is not in the
// vendored subset)
struct InputSegment {
  unsigned int firstLine;
  unsigned int firstFileLine;
  std::string fileName;
};


struct PSMDescriptionPriority {
  PSMDescription* psm;
  size_t priority;
  LabelType label;
  
  bool operator<(const PSMDescriptionPriority& psmp) const {
    return (priority < psmp.priority);
  }
};

/*
* ScanId is a tuple of (ScanNr, ExpMass).
*/
typedef std::pair<int, double> ScanId;

/*
* What reading the input found, besides the PSMs, for the caller to pass on
* (to the SanityCheck of the run).
*/
struct InputSummary {
  // No spectrum (scan number and experimental mass) has both a target and a
  // decoy PSM: the input comes from a concatenated search rather than from
  // separate target and decoy searches.
  bool concatenatedSearch = true;
  // the initial direction of the DefaultDirection line; empty if the input
  // has none or if it is all 0
  std::vector<double> defaultWeights;
};

/*
* SetHandler is a class that provides functionality to handle training,
* testing, Xval data sets, reads/writes from/to a file, prints them.
*
*/
class SetHandler {
 public:
  SetHandler(unsigned int maxPSMs);
  virtual ~SetHandler();
  SetHandler(const SetHandler&) = delete;
  SetHandler& operator=(const SetHandler&) = delete;

  void push_back_dataset(std::unique_ptr<DataSet> ds);

  void setDecoyPrefix(std::string& decoyPrefix) { decoyPrefix_ = decoyPrefix; }

     
  //const double* getFeatures(const int setPos, const int ixPos) const; 
  size_t getMaxPSMs() { return maxPSMs_; }
  // number of PSM lines in the input, counted by the last full read
  size_t getNumPsmLines() const { return numPsmLines_; }
  
  // Reads the PSMs of a tab delimited stream into the target and the decoy
  // set and describes the input in summary. Returns 0 on error, without
  // changing summary, and 1 on success.
  int readTab(std::istream& dataStream, InputSummary& summary);
  // The second pass of -N: scores every PSM of the stream with rawWeights
  // and adds it to allScores; a PSM whose input line is listed in lineFolds
  // (sorted by line) is scored with foldWeights[fold] instead.
  int readAndScoreTab(std::istream& dataStream,
    std::vector<double>& rawWeights, Scores& allScores,
    const std::vector<std::vector<double> >& foldWeights = {},
    const std::vector<std::pair<unsigned int, unsigned int> >& lineFolds = {});
  void addQueueToSets(std::priority_queue<PSMDescriptionPriority>& subsetPSMs,
    DataSet* targetSet, DataSet* decoySet);
  
  // defaultWeights: the initial direction to write as the DefaultDirection
  // line, if not empty
  void writeTab(const std::string& dataFN,
                const std::vector<double>& defaultWeights);
  void populateScoresWithPSMs(std::vector<ScoreHolder> &scores, LabelType label);
  void normalizeFeatures(Normalizer& norm);
  
  LabelType getLabel(int setPos);
  inline int getSizeFromLabel(LabelType label) {
    return static_cast<int>(subsets_[getSubsetIndexFromLabel(label)]->getSize());
  }
  
  inline DataSet* getSubset(unsigned int ix) { return subsets_[ix].get(); }
  inline DataSet* getSubsetFromLabel(LabelType label) {
    return subsets_[getSubsetIndexFromLabel(label)].get();
  }
  
  FeatureMemoryPool& getFeaturePool() { return featurePool_; }

  // For an input concatenated from several files: where the lines of each
  // file are, so that messages can name the file and line of a PSM.
  void setInputSegments(std::vector<InputSegment> segments) {
    ownSegments_ = std::move(segments);
    inputSegments_ = &ownSegments_;
  }
  // The same for a list that the reader of the files extends while the
  // lines are read (ConcatenatedPins::segments); it has to outlive the
  // reading.
  void setInputSegments(const std::vector<InputSegment>* segments) {
    inputSegments_ = segments;
  }
  // "line <lineNr>", or "line <n> of <file>" for a concatenated input
  std::string describeLine(unsigned int lineNr) const;

  // Among the spectra (ScanNr and ExpMass) of the first kRepeatSampleLines
  // PSM lines of the last read: their number, and the number with more
  // than one PSM of the same label, as in a search that reports several
  // candidates per spectrum.
  size_t getSampledLines() const { return sampledLines_; }
  size_t getSampledSpectra() const { return sampledSpectra_; }
  size_t getSampledSpectraWithRepeatedLabel() const {
    return sampledSpectraWithRepeatedLabel_;
  }
  static constexpr size_t kRepeatSampleLines = 100000u;
  
  void reset();
  void releaseFeatures();

 protected:
  size_t maxPSMs_;
  size_t numPsmLines_ = 0;
  PsmStore psmStore_;  // owns the PSMs of subsets_ and of -N pass 2
  std::vector<std::unique_ptr<DataSet> > subsets_;
  FeatureMemoryPool featurePool_;
  std::string decoyPrefix_; // Used to determine if a psm is a decoy
  // empty for a single input file
  std::vector<InputSegment> ownSegments_;
  const std::vector<InputSegment>* inputSegments_ = &ownSegments_;
  size_t sampledLines_ = 0u, sampledSpectra_ = 0u,
         sampledSpectraWithRepeatedLabel_ = 0u;
  // labels seen per spectrum while reading: target (any label but -1),
  // decoy, and whether one of them was seen twice
  static constexpr unsigned char kSeenTarget = 1u, kSeenDecoy = 2u,
                                 kSeenTwice = 4u;
  static unsigned char labelBit(bool isDecoy) {
    return isDecoy ? kSeenDecoy : kSeenTarget;
  }
  // adds a label to those seen for a spectrum, counts a repeated label
  // while sampling, and returns false if the spectrum now has both labels
  bool addLabel(unsigned char& seen, bool isDecoy, bool sampling);
  // warnings about labels other than 1 and -1 printed per line, then counted
  static constexpr unsigned int kMaxLabelWarnings = 1000u;
  
  // rethrows e, with the file and line of a PSM line of a concatenated input
  [[noreturn]] void rethrowWithSource(const MyException& e,
                                      unsigned int lineNr) const;
  // Moves on to the next input line. Line numbers, and with them the
  // default scan numbers and the PSM counts, are kept in 32-bit integers,
  // also signed ones, so an input with more than INT_MAX lines is an error.
  static void incrementLineNr(unsigned int& lineNr);
  void warnIgnoredLabel(unsigned int lineNr, unsigned int& numIgnored) const;
  void warnIgnoredLabelTotal(unsigned int numIgnored) const;
  
  // readTab (summary given, no weights) and readAndScoreTab
  int readOrScoreTab(std::istream& dataStream,
    std::vector<double>& rawWeights, Scores& allScores,
    const std::vector<std::vector<double> >& foldWeights,
    const std::vector<std::pair<unsigned int, unsigned int> >& lineFolds,
    InputSummary* summary);

  unsigned int getSubsetIndexFromLabel(LabelType label);
  static inline std::string &rtrim(std::string &s);
  
  int getOptionalFields(const std::string& headerLine, 
    std::vector<OptionalField>& optionalFields);
  bool isDefaultDirectionLine(const std::string& defaultDirectionLine);
  int getNumFeatures(const std::string& line, int optionalFieldCount);
  void getFeatureNames(const std::string& headerLine, int numFeatures, 
    int optionalFieldCount, FeatureNames& featureNames);
  void checkPeptideColumn(const std::string& headerLine, int numFeatures,
    int optionalFieldCount, unsigned int firstPsmLineNr) const;
  bool getInitValues(const std::string& defaultDirectionLine, 
    int optionalFieldCount, std::vector<double>& init_values);
  ScanId getScanId(const std::string& psmLine, int& label,
    std::vector<OptionalField>& optionalFields, unsigned int lineNr);
  bool readScanIdAndLabel(const std::string& psmLine, unsigned int lineNr,
    std::vector<OptionalField>& optionalFields, TabReader& reader,
    PsmHeader& header, ScanId& scanId, int& label);
    
  // the loop over the PSM lines and the reader of one line (SetHandler.cpp)
  template <class PerLine>
  unsigned int forEachPsmLine(std::istream& dataStream, std::string& psmLine,
    unsigned int firstLineNr, const char* progressVerb, PerLine perLine);
  class PsmLineReader;

  void readPSMs(std::istream& dataStream, std::string& psmLine,
    bool hasInitialValueRow, bool& separateSearches,
    std::vector<OptionalField>& optionalFields);
  static void rescaleFoldScores(Scores& allScores,
    const std::vector<std::pair<size_t, unsigned int> >& rows, size_t numFolds);
  void readAndScorePSMs(std::istream& dataStream, std::string& psmLine,
    bool hasInitialValueRow, std::vector<OptionalField>& optionalFields,
    std::vector<double>& rawWeights, Scores& allScores,
    const std::vector<std::vector<double> >& foldWeights,
    const std::vector<std::pair<unsigned int, unsigned int> >& lineFolds);
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*SETHANDLER_H_*/
