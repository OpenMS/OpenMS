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
    @brief Global-alignment similarity score for two amino-acid sequences using the Needleman-Wunsch algorithm.

    Computes the score (not the alignment itself) of the best global
    alignment of two amino-acid sequences against a substitution matrix
    plus a linear gap penalty. The implementation uses a two-row
    rolling buffer.

    @note The implementation indexes the substitution matrix by
          @c c - 'A', so the input strings must contain only uppercase
          ASCII letters in @c [A-Z]. Characters @c J, @c O, and @c U
          have @c INT16_MAX rows/columns in both bundled matrices and
          effectively forbid any alignment involving them. Lowercase or
          non-ASCII characters lead to out-of-bounds reads.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI NeedlemanWunsch
  {

  public:

    /// Substitution matrices bundled with this class.
    enum class ScoringMatrix
    {
      identity,                ///< +1 on the diagonal, 0 elsewhere; @c J / @c O / @c U columns and rows are forbidden (@c INT16_MAX).
      PAM30MS,                 ///< PAM30 matrix extended for MS-style substitutions (default).
      SIZE_OF_SCORINGMATRIX    ///< Number of valid matrices.
    };

    /**
      @brief Construct with a chosen scoring matrix and gap penalty.

      @param[in] matrix  Scoring matrix to use.
      @param[in] penalty Linear gap penalty (subtracted per gap step).
    */
    NeedlemanWunsch(ScoringMatrix matrix, int penalty);

    /**
      @brief Default constructor.

      Uses @c ScoringMatrix::PAM30MS and a gap penalty of @c 5.
    */
    NeedlemanWunsch() = default;

    /**
      @brief Destructor.
    */
    ~NeedlemanWunsch()=default;

    /// Lookup table for @ref setMatrix(const std::string&); current entries: @c "identity", @c "PAM30MS".
    static const std::vector<std::string> NamesOfScoringMatrices;

    /**
      @brief Compute the Needleman-Wunsch similarity score of @p seq1 vs. @p seq2.

      Uses the currently selected scoring matrix and gap penalty. The
      sequences must contain only uppercase @c [A-Z] characters; see
      the class brief for the constraints.

      @param[in] seq1 First amino-acid sequence.
      @param[in] seq2 Second amino-acid sequence.
      @return Score of the best global alignment.
    */
    int align(const String& seq1, const String& seq2);

    /**
      @brief Select the scoring matrix used by @ref align.

      @param[in] matrix Scoring matrix to use.
    */
    void setMatrix(const ScoringMatrix& matrix);

    /**
      @brief Select the scoring matrix by name; the accepted names are listed in @ref NamesOfScoringMatrices.

      @param[in] matrix Matrix name (currently @c "identity" or @c "PAM30MS").

      @throws Exception::IllegalArgument when @p matrix is not one of
                                         the names in
                                         @ref NamesOfScoringMatrices.
                                         The error message lists the
                                         valid names.
    */
    void setMatrix(const std::string& matrix);

    /**
      @brief Set the linear gap penalty subtracted per gap step.

      @param[in] penalty Gap penalty (typically positive); the
                         algorithm subtracts @p penalty from the score
                         for each gap step.
    */
    void setPenalty(const int penalty);

    /**
      @brief Currently selected scoring matrix.

      @return The scoring matrix passed to the last @ref setMatrix call
              (or the default @c PAM30MS if @ref setMatrix has not been
              called).
    */
    ScoringMatrix getMatrix() const;

    /**
      @brief Currently configured gap penalty.

      @return The penalty passed to the last @ref setPenalty call (or
              the default @c 5 if @ref setPenalty has not been called).
    */
    int getPenalty() const;

  private:
    int gap_penalty_ = 5; ///< Linear gap penalty.
    ScoringMatrix my_matrix_ = ScoringMatrix::PAM30MS; ///< Currently selected scoring matrix.
    std::vector<int> first_row_{};  ///< First row of the two-row rolling buffer used by @ref align.
    std::vector<int> second_row_{}; ///< Second row of the two-row rolling buffer used by @ref align.
  };

}