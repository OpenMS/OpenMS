// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

    defaults_.setValue("matrix", "PAM30MS", "Substitution matrix to use for alignment-based similarity scoring");
    defaults_.setValidStrings("matrix", NeedlemanWunsch::NamesOfScoringMatrices);
    defaults_.setValue("penalty", 5, "Alignment gap penalty (the same value is used for gap opening and extension)");
    defaults_.setMinInt("penalty", 1);

    defaultsToParam_();

  }

  void ConsensusIDAlgorithmPEPMatrix::updateMembers_()
  {
    ConsensusIDAlgorithmSimilarity::updateMembers_();

    string matrix = param_.getValue("matrix");
    int penalty = param_.getValue("penalty");

    alignment_.setMatrix(matrix);

    if (penalty > 0)
    {
      alignment_.setPenalty(penalty);
    }
    else
    {
      String msg = "Gap penalty should be positive";
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

    double score_sim = alignment_.align(unmod_seq1, unmod_seq2);

    if (score_sim < 0)
    {
      score_sim = 0;
    }
    else
    {
      double score_self1 = alignment_.align(unmod_seq1, unmod_seq1);
      double score_self2 = alignment_.align(unmod_seq2, unmod_seq2);
      score_sim /= min(score_self1, score_self2); // normalize
    }
    return score_sim;
  }

} // namespace OpenMS
