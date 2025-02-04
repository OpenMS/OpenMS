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
    analysis_results_(nullptr),
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
      analysis_results_(nullptr),
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
    analysis_results_(nullptr),
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
    analysis_results_(nullptr),
    rank_(source.rank_),
    charge_(source.charge_),
    peptide_evidences_(source.peptide_evidences_),
    fragment_annotations_(source.fragment_annotations_)
  {
    if (source.analysis_results_ != nullptr)
    {
      analysis_results_ = new std::vector<PepXMLAnalysisResult>(*source.analysis_results_);
    }
  }

  /// Move constructor
  PeptideHit::PeptideHit(PeptideHit&& source) noexcept :
    MetaInfoInterface(std::move(source)), // NOTE: rhs itself is an lvalue
    sequence_(std::move(source.sequence_)),
    score_(source.score_),
    analysis_results_(std::move(source.analysis_results_)),
    rank_(source.rank_),
    charge_(source.charge_),
    peptide_evidences_(std::move(source.peptide_evidences_)),
    fragment_annotations_(std::move(source.fragment_annotations_))
  {
    // see http://thbecker.net/articles/rvalue_references/section_05.html
    source.analysis_results_ = nullptr;
  }

  // destructor
  PeptideHit::~PeptideHit()
  {
    delete analysis_results_;
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
    delete analysis_results_;
    if (source.analysis_results_ != nullptr)
    {
      analysis_results_ = new std::vector<PepXMLAnalysisResult>(*source.analysis_results_);
    }
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

    // free memory and assign rhs memory
    delete analysis_results_;
    analysis_results_ = source.analysis_results_;
    source.analysis_results_ = nullptr;

    rank_ = source.rank_;
    charge_ = source.charge_;
    peptide_evidences_ = source.peptide_evidences_;
    fragment_annotations_ = source.fragment_annotations_;

    return *this;
  }

  bool PeptideHit::operator==(const PeptideHit& rhs) const
  {
    bool ar_equal = false;
    if (analysis_results_ == nullptr && rhs.analysis_results_ == nullptr)
    {
      ar_equal = true;
    }
    else if (analysis_results_ != nullptr && rhs.analysis_results_ != nullptr)
    {
      ar_equal = (*analysis_results_ == *rhs.analysis_results_);
    }
    else
    {
      return false; // one is null the other isn't
    }
    return MetaInfoInterface::operator==(rhs)
           && sequence_ == rhs.sequence_
           && score_ == rhs.score_
           && ar_equal
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

  void PeptideHit::setAnalysisResults(std::vector<PeptideHit::PepXMLAnalysisResult> aresult)
  {
    // delete old results first
    if (analysis_results_ != nullptr) delete analysis_results_;
    analysis_results_ = new std::vector< PeptideHit::PepXMLAnalysisResult> (std::move(aresult));
  }

  void PeptideHit::addAnalysisResults(const PeptideHit::PepXMLAnalysisResult& aresult)
  {
    if (analysis_results_ == nullptr)
    {
      analysis_results_ = new std::vector< PeptideHit::PepXMLAnalysisResult>();
    }
    analysis_results_->push_back(aresult);
  }
  
  const std::vector<PeptideHit::PepXMLAnalysisResult>& PeptideHit::getAnalysisResults() const
  {
    static std::vector<PeptideHit::PepXMLAnalysisResult> empty;
    if (analysis_results_ == nullptr)
    {
      return empty;
    }
    return (*analysis_results_);
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
