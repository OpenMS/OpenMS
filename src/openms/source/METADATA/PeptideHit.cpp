// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/PeptideHit.h>

#include <algorithm>

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

  // destructor
  PeptideHit::~PeptideHit()
  {
    if (analysis_results_ != nullptr)
    {
      // free memory again
      delete analysis_results_;
    }
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
    analysis_results_ = nullptr;
    if (source.analysis_results_ != nullptr)
    {
      if (analysis_results_ != nullptr)
      {
        // free memory first
        delete analysis_results_;
      }
      analysis_results_ = new std::vector<PepXMLAnalysisResult>(*source.analysis_results_);
    }
    charge_ = source.charge_;
    rank_  = source.rank_;
    peptide_evidences_ = source.peptide_evidences_;
    fragment_annotations_ = source.fragment_annotations_;
    return *this;
  }

  bool PeptideHit::operator==(const PeptideHit& rhs) const
  {
    bool ar_equal = false;
    if (analysis_results_ == nullptr && rhs.analysis_results_ == nullptr) ar_equal = true;
    else if (analysis_results_ != nullptr && rhs.analysis_results_ != nullptr)
    {
      ar_equal = (*analysis_results_ == *rhs.analysis_results_);
    }
    else return false; // one is null the other isn't

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

  // returns the peptide sequence without trailing or following spaces
  const AASequence& PeptideHit::getSequence() const
  {
    return sequence_;
  }

  void PeptideHit::setSequence(const AASequence& sequence)
  {
    sequence_ = sequence;
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
    analysis_results_ = new std::vector< PeptideHit::PepXMLAnalysisResult> (aresult);
  }

  void PeptideHit::addAnalysisResults(PeptideHit::PepXMLAnalysisResult aresult)
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
    for (vector<PeptideEvidence>::const_iterator it = peptide_evidences_.begin(); it != peptide_evidences_.end(); ++it)
    {
      // don't return empty accessions
      if (!it->getProteinAccession().empty())
      {
        accessions.insert(it->getProteinAccession());
      }
    }
    return accessions;
  }

  std::vector<PeptideHit::PeakAnnotation> PeptideHit::getPeakAnnotations() const
  {
    return fragment_annotations_;
  }

  void PeptideHit::setPeakAnnotations(std::vector<PeptideHit::PeakAnnotation> frag_annotations)
  {
    fragment_annotations_ = frag_annotations;
  }

} // namespace OpenMS
