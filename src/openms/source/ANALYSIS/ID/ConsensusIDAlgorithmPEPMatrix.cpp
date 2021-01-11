// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmPEPMatrix.h>

using namespace std;

namespace OpenMS
{
  ConsensusIDAlgorithmPEPMatrix::ConsensusIDAlgorithmPEPMatrix()
  {
    setName("ConsensusIDAlgorithmPEPMatrix"); // DefaultParamHandler

    defaults_.setValue("matrix", "identity", "Substitution matrix to use for alignment-based similarity scoring");
    defaults_.setValidStrings("matrix", ListUtils::create<String>("identity,PAM30MS"));
    defaults_.setValue("penalty", 5, "Alignment gap penalty (the same value is used for gap opening and extension)");
    defaults_.setMinInt("penalty", 1);

    defaultsToParam_();

    ::seqan::resize(rows(alignment_), 2);
  }


  void ConsensusIDAlgorithmPEPMatrix::updateMembers_()
  {
    ConsensusIDAlgorithmSimilarity::updateMembers_();

    // alignment scoring using SeqAn/similarity matrices:
    String matrix = param_.getValue("matrix");
    int penalty = param_.getValue("penalty");
    scoring_method_ = SeqAnScore(-penalty, -penalty);
    if (matrix == "identity")
    {
      ::seqan::setDefaultScoreMatrix(scoring_method_, 
                                     ::seqan::AdaptedIdentity());
    }
    else if (matrix == "PAM30MS")
    {
      ::seqan::setDefaultScoreMatrix(scoring_method_, ::seqan::PAM30MS());
    }
    else
    {
      String msg = "Matrix '" + matrix + "' is not known! Valid choices are: "
        "'identity', 'PAM30MS'.";
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       msg);
    }

    // new parameters may affect the similarity calculation, so clear cache:
    similarities_.clear();
  }


  double ConsensusIDAlgorithmPEPMatrix::getSimilarity_(AASequence seq1,
                                                       AASequence seq2)
  {
    // here we cannot take modifications into account:
    String unmod_seq1 = seq1.toUnmodifiedString();
    String unmod_seq2 = seq2.toUnmodifiedString();
    if (unmod_seq1 == unmod_seq2) return 1.0;
    // order of sequences matters for cache look-up:
    if (unmod_seq1 > unmod_seq2) swap(unmod_seq1, unmod_seq2);
    seq1 = AASequence::fromString(unmod_seq1);
    seq2 = AASequence::fromString(unmod_seq2);
    pair<AASequence, AASequence> seq_pair = make_pair(seq1, seq2);
    SimilarityCache::iterator pos = similarities_.find(seq_pair);
    if (pos != similarities_.end()) return pos->second; // score found in cache
    
    // use SeqAn similarity scoring:
    SeqAnSequence seqan_seq1 = unmod_seq1.c_str();
    SeqAnSequence seqan_seq2 = unmod_seq2.c_str();
    // seq. 1 against itself:
    ::seqan::assignSource(row(alignment_, 0), seqan_seq1);
    ::seqan::assignSource(row(alignment_, 1), seqan_seq1);
    double score_self1 = globalAlignment(alignment_, scoring_method_,
                                         ::seqan::NeedlemanWunsch());
    // seq. 1 against seq. 2:
    ::seqan::assignSource(row(alignment_, 1), seqan_seq2);
    double score_sim = globalAlignment(alignment_, scoring_method_, 
                                       ::seqan::NeedlemanWunsch());
    // seq. 2 against itself:
    ::seqan::assignSource(row(alignment_, 0), seqan_seq2);
    double score_self2 = globalAlignment(alignment_, scoring_method_,
                                         ::seqan::NeedlemanWunsch());
    if (score_sim < 0)
    {
      score_sim = 0;
    }
    else
    {
      score_sim /= min(score_self1, score_self2); // normalize
    }
    similarities_[seq_pair] = score_sim; // cache the similarity score

    return score_sim;
  }

} // namespace OpenMS
