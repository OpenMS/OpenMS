// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideHit.h>
#include <ostream>
#include <utility>

using namespace std;

namespace OpenMS
{
  // default constructor
  PeptideHit::PeptideHit() :
    MetaInfoInterface(),
    sequence_(),
    score_(0),
    rank_(0),
    charge_(0),
    peptide_evidences_(),
    fragment_annotations_()
  {
  }

  // values constructor
  PeptideHit::PeptideHit(double score, UInt rank, Int charge, const AASequence& sequence) :
      MetaInfoInterface(),
      sequence_(sequence),
      score_(score),
      rank_(rank),
      charge_(charge),
      peptide_evidences_(),
      fragment_annotations_()
  {
  }

  // values constructor
  PeptideHit::PeptideHit(double score, UInt rank, Int charge, AASequence&& sequence) :
    MetaInfoInterface(),
    sequence_(sequence),
    score_(score),
    rank_(rank),
    charge_(charge),
    peptide_evidences_(),
    fragment_annotations_()
  {
  }

  // copy constructor
  PeptideHit::PeptideHit(const PeptideHit& source) :
    MetaInfoInterface(source),
    sequence_(source.sequence_),
    score_(source.score_),
    rank_(source.rank_),
    charge_(source.charge_),
    peptide_evidences_(source.peptide_evidences_),
    fragment_annotations_(source.fragment_annotations_)
  {
  }

  /// Move constructor
  PeptideHit::PeptideHit(PeptideHit&& source) noexcept :
    MetaInfoInterface(std::move(source)), // NOTE: rhs itself is an lvalue
    sequence_(std::move(source.sequence_)),
    score_(source.score_),
    rank_(source.rank_),
    charge_(source.charge_),
    peptide_evidences_(std::move(source.peptide_evidences_)),
    fragment_annotations_(std::move(source.fragment_annotations_))
  {
  }

  // destructor
  PeptideHit::~PeptideHit()
  {
  }

  PeptideHit& PeptideHit::operator=(const PeptideHit& source)
  {
    if (this == &source)
    {
      return *this;
    }

    MetaInfoInterface::operator=(source);
    sequence_ = source.sequence_;
    score_ = source.score_;
    rank_ = source.rank_;
    charge_ = source.charge_;
    peptide_evidences_ = source.peptide_evidences_;
    fragment_annotations_ = source.fragment_annotations_;
    return *this;
  }

  PeptideHit& PeptideHit::operator=(PeptideHit&& source) noexcept
  {
    if (&source == this)
    {
      return *this;
    }

    MetaInfoInterface::operator=(std::move(source));
    //clang-tidy overly strict, should be fine to move the rest here
    sequence_ = source.sequence_;
    score_ = source.score_;
    rank_ = source.rank_;
    charge_ = source.charge_;
    peptide_evidences_ = source.peptide_evidences_;
    fragment_annotations_ = source.fragment_annotations_;

    return *this;
  }

  bool PeptideHit::operator==(const PeptideHit& rhs) const
  {
    return MetaInfoInterface::operator==(rhs)
           && sequence_ == rhs.sequence_
           && score_ == rhs.score_
           && rank_ == rhs.rank_
           && charge_ == rhs.charge_
           && peptide_evidences_ == rhs.peptide_evidences_
           && fragment_annotations_ == rhs.fragment_annotations_;
  }

  bool PeptideHit::operator!=(const PeptideHit& rhs) const
  {
    return !operator==(rhs);
  }

  // returns the score of the peptide hit
  double PeptideHit::getScore() const
  {
    return score_;
  }

  // returns the rank of the peptide hit
  UInt PeptideHit::getRank() const
  {
    return rank_;
  }

  const AASequence& PeptideHit::getSequence() const
  {
    return sequence_;
  }

  AASequence& PeptideHit::getSequence()
  {
    return sequence_;
  }

  void PeptideHit::setSequence(const AASequence& sequence)
  {
    sequence_ = sequence;
  }

  void PeptideHit::setSequence(AASequence&& sequence)
  {
    sequence_ = std::move(sequence);
  }

  Int PeptideHit::getCharge() const
  {
    return charge_;
  }

  void PeptideHit::setCharge(Int charge)
  {
    charge_ = charge;
  }

  const std::vector<PeptideEvidence>& PeptideHit::getPeptideEvidences() const
  {
    return peptide_evidences_;
  }

  void PeptideHit::setPeptideEvidences(const std::vector<PeptideEvidence>& peptide_evidences)
  {
    peptide_evidences_ = peptide_evidences;
  }

  void PeptideHit::setPeptideEvidences(std::vector<PeptideEvidence>&& peptide_evidences)
  {
    peptide_evidences_ = std::move(peptide_evidences);
  }

  void PeptideHit::addPeptideEvidence(const PeptideEvidence& peptide_evidence)
  {
    peptide_evidences_.push_back(peptide_evidence);
  }

  // sets the score of the peptide hit
  void PeptideHit::setScore(double score)
  {
    score_ = score;
  }

  void PeptideHit::setAnalysisResults(const std::vector<PeptideHit::PepXMLAnalysisResult>& aresult)
  {
    // Remove all existing analysis result meta values
    std::vector<String> keys;
    getKeys(keys);
    for (const auto& key : keys)
    {
      if (key.hasPrefix("_ar_"))
      {
        removeMetaValue(key);
      }
    }
    
    // Add new analysis results as meta values
    for (size_t i = 0; i < aresult.size(); ++i)
    {
      const auto& ar = aresult[i];
      setMetaValue("_ar_" + String(i) + "_score_type", ar.score_type);
      setMetaValue("_ar_" + String(i) + "_score", ar.main_score);
      setMetaValue("_ar_" + String(i) + "_higher_is_better", ar.higher_is_better == true ? "true" : "false");
      
      for (const auto& subscore : ar.sub_scores)
      {
        setMetaValue("_ar_" + String(i) + "_subscore_" + subscore.first, subscore.second);
      }
    }
  }

  void PeptideHit::addAnalysisResults(const PeptideHit::PepXMLAnalysisResult& aresult)
  {
    size_t index = getNumberOfAnalysisResultsFromMetaValues_();
    
    setMetaValue("_ar_" + String(index) + "_score_type", aresult.score_type);
    setMetaValue("_ar_" + String(index) + "_score", aresult.main_score);
    setMetaValue("_ar_" + String(index) + "_higher_is_better", aresult.higher_is_better == true ? "true" : "false");
    
    for (const auto& subscore : aresult.sub_scores)
    {
      setMetaValue("_ar_" + String(index) + "_subscore_" + subscore.first, subscore.second);
    }
  }
  
  std::vector<PeptideHit::PepXMLAnalysisResult> PeptideHit::getAnalysisResults() const
  {
    return extractAnalysisResultsFromMetaValues_();
  }

  size_t PeptideHit::getNumberOfAnalysisResultsFromMetaValues_() const
  {
    size_t count = 0;
    std::vector<String> keys;
    getKeys(keys);
    
    for (const auto& key : keys)
    {
      if (key.hasPrefix("_ar_") &&
          key.hasSuffix("_score_type"))
      {
        ++count;
      }
    }
    
    return count;
  }

  std::vector<PeptideHit::PepXMLAnalysisResult> PeptideHit::extractAnalysisResultsFromMetaValues_() const
  {
    std::vector<PeptideHit::PepXMLAnalysisResult> results;
    std::vector<String> keys;
    getKeys(keys);
    
    // First, find all indices that have analysis results
    std::set<size_t> indices;
    
    for (const auto& key : keys)
    {
      const String prefix = "_ar_";
      const String suffix = "_score_type";
      if (key.hasPrefix(prefix) &&
          key.hasSuffix(suffix))
      {
        String index_str = key.substr(prefix.size(), key.size() - prefix.size() - suffix.size()); // Extract index from _ar_<index>_score_type"
        indices.insert(index_str.toInt());
      }
    }
    
    // For each index, extract the analysis result
    for (size_t index : indices)
    {
      PeptideHit::PepXMLAnalysisResult ar;
      String prefix = "_ar_" + String(index) + "_";
      
      // Get score type
      if (metaValueExists(prefix + "score_type"))
      {
        ar.score_type = getMetaValue(prefix + "score_type").toString();
      }
      
      // Get main score
      if (metaValueExists(prefix + "score"))
      {
        ar.main_score = getMetaValue(prefix + "score");
      }
      
      // Get higher_is_better flag
      if (metaValueExists(prefix + "higher_is_better"))
      {
        ar.higher_is_better = getMetaValue(prefix + "higher_is_better").toBool();
      }
      
      // Get sub-scores
      String subscore_prefix = prefix + "subscore_";
      for (const auto& key : keys)
      {
        if (key.hasPrefix(subscore_prefix))
        {
          String subscore_name = key.substr(subscore_prefix.size());
          ar.sub_scores[subscore_name] = getMetaValue(key);
        }
      }
      
      results.push_back(ar);
    }
    
    return results;
  }

  // sets the rank
  void PeptideHit::setRank(UInt newrank)
  {
    rank_ = newrank;
  }

  std::set<String> PeptideHit::extractProteinAccessionsSet() const
  {
    set<String> accessions;
    for (const auto& ev : peptide_evidences_)
    {
      // don't return empty accessions
      if (!ev.getProteinAccession().empty())
      {
        accessions.insert(ev.getProteinAccession());
      }
    }
    return accessions;
  }

  std::vector<PeptideHit::PeakAnnotation>& PeptideHit::getPeakAnnotations()
  {
    return fragment_annotations_;
  }

  const std::vector<PeptideHit::PeakAnnotation>& PeptideHit::getPeakAnnotations() const
  {
    return fragment_annotations_;
  }

  void PeptideHit::setPeakAnnotations(std::vector<PeptideHit::PeakAnnotation> frag_annotations)
  {
    fragment_annotations_ = std::move(frag_annotations);
  }

  std::ostream& operator<< (std::ostream& stream, const PeptideHit& hit)
  {
    return stream << "peptide hit with sequence '" + hit.getSequence().toString() +
           "', charge " + String(hit.getCharge()) + ", score " +
           String(hit.getScore());
  }

} // namespace OpenMS
