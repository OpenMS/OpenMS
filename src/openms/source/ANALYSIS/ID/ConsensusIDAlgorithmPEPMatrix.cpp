// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
