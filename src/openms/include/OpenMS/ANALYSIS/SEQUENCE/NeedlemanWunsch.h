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
// $Maintainer: Nora Wild $
// $Authors: Nora Wild $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <vector>

namespace OpenMS
{
  /**
   @brief This class serves for the calculation of the global alignment score of two amino acid sequences
   by using the Needleman-Wunsch-Algorithm. For match and mismatch it uses a similarity scoring matrix.
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