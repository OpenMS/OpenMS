// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Nora Wild $
// $Authors: Nora Wild $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <vector>

namespace OpenMS
{
  /**
   @brief This class contains functions that are used to calculate the global alignment score of two amino acid sequences. 
   This class uses the Needleman-Wunsch algorithm. For match and mismatch it uses a similarity scoring matrix.
   */
  class OPENMS_DLLAPI NeedlemanWunsch
  {

  public:

    /// contains the valid matrices and the number of them
    enum class ScoringMatrix
    {
      identity,
      PAM30MS,
      SIZE_OF_SCORINGMATRIX
    };

    /// Constructor that sets the scoring matrix and the gap penalty
    NeedlemanWunsch(ScoringMatrix matrix, int penalty);

    /// Default constructor (scoring matrix PAM30MS and penalty 5)
    NeedlemanWunsch() = default;

    /// Default destructor
    ~NeedlemanWunsch()=default;

    /// Names of valid matrices
    static const std::vector<std::string> NamesOfScoringMatrices;

    /**
     @brief Calculates the similarity score of the global alignment of two amino acid sequences using Needleman-Wunsch-
     Algorithm.
     */
    int align(const String& seq1, const String& seq2);

    /**
     @brief sets the scoring matrix. Takes either a string or the enum ScoringMatrix.
     @exception Exception: illegal argument is thrown if the input is not a member of the valid matrices.
     */
    void setMatrix(const ScoringMatrix& matrix);
    void setMatrix(const std::string& matrix);

    /// sets the cost of gaps
    void setPenalty(const int penalty);

    /// returns the scoring matrix
    ScoringMatrix getMatrix() const;

    /// returns the gap penalty
    int getPenalty() const;

  private:
    int gap_penalty_ = 5; ///< penalty for alignment score calculation
    ScoringMatrix my_matrix_ = ScoringMatrix::PAM30MS; ///< scoring matrix for the alignment score calculation
    std::vector<int> first_row_{}; ///< alignment score calculation with two rows
    std::vector<int> second_row_{};
  };

}