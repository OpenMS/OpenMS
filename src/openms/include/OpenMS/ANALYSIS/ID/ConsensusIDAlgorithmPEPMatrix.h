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
// $Maintainer: Hendrik Weisser $
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMPEPMATRIX_H
#define OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMPEPMATRIX_H

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>

// Extend SeqAn by a user-define scoring matrix.
namespace seqan
{

  // We have to create a new specialization of the _ScoringMatrix class
  // for amino acids.  For this, we first create a new tag.
  struct PAM30MS {}; // PAM30MS matrix
  struct AdaptedIdentity {}; // identity matrix adapted for I/L, Q/K ambiguity

  // Then, we specialize the class _ScoringMatrix.
  template <>
  struct ScoringMatrixData_<int, AminoAcid, PAM30MS>
  {
    enum
    {
      VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
      TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };
    static inline const int* getData()
    {
      // Rant: I cannot find a primary source for the PAM30MS scoring matrix!
      // It seems to have been first published in Huang et al., JBC 2001
      // (http://www.jbc.org/content/276/30/28327), but the paper does not show
      // the actual matrix (gah!).
      // The matrix here comes from old OpenMS code and also matches this one:
      // http://proteomics.fiocruz.br/supplementaryfiles/pepexplorer/BeforeRevision/PFUGridResults/PFUGridSearch/pam30ms.txt

      static const int _data[TAB_SIZE] =
      {
        //        A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        /* A */   6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2,  0, -1,-13, -8, -2, -7, -6,  0,-17,
        /* R */  -7,  8, -6,-10, -8, -2, -9, -9, -2, -5, -7,  0, -4, -9, -4, -3, -6, -2,-10, -8,  5, -1,  0,-17,
        /* N */  -4, -6,  8,  2,-11, -3, -2, -3,  0, -5, -6, -1, -9, -9, -6,  0, -2, -8, -4, -8, -4, -2,  0,-17,
        /* D */  -3,-10,  2,  8,-14, -2,  2, -3, -4, -7,-10, -4,-11,-15, -8, -4, -5,-15,-11, -8, -7, -3,  0,-17,
        /* C */  -6, -8,-11,-14, 10,-14,-14, -9, -7, -6,-11,-14,-13,-13, -8, -3, -8,-15, -4, -6,-11,-14,  0,-17,
        /* Q */  -4, -2, -3, -2,-14,  8,  1, -7,  1, -8, -7, -3, -4,-13, -3, -5, -5,-13,-12, -7, -3,  4,  0,-17,
        /* E */  -2, -9, -2,  2,-14,  1,  8, -4, -5, -5, -7, -4, -7,-14, -5, -4, -6,-17, -8, -6, -7, -2,  0,-17,
        /* G */  -2, -9, -3, -3, -9, -7, -4,  6, -9,-11,-11, -7, -8, -9, -6, -2, -6,-15,-14, -5, -8, -7,  0,-17,
        /* H */  -7, -2,  0, -4, -7,  1, -5, -9,  9, -9, -8, -6,-10, -6, -4, -6, -7, -7, -3, -6, -4, -3,  0,-17,
        /* I */  -5, -5, -5, -7, -6, -8, -5,-11, -9,  8,  5, -6, -1, -2, -8, -7, -2,-14, -6,  2, -6, -7,  0,-17,
        /* L */  -6, -7, -6,-10,-11, -7, -7,-11, -8,  5,  5, -7,  0, -3, -8, -8, -5,-10, -7,  0, -7, -7,  0,-17,
        /* K */  -7,  0, -1, -4,-14, -3, -4, -7, -6, -6, -7,  7, -2,-14, -6, -4, -3,-12, -9, -9,  5,  4,  0,-17,
        /* M */  -5, -4, -9,-11,-13, -4, -7, -8,-10, -1,  0, -2, 11, -4, -8, -5, -4,-13,-11, -1, -3, -3,  0,-17,
        /* F */  -8, -9, -9,-15,-13,-13,-14, -9, -6, -2, -3,-14, -4,  9,-10, -6, -9, -4,  2, -8,-12,-14,  0,-17,
        /* P */  -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -8, -6, -8,-10,  8, -2, -4,-14,-13, -6, -5, -5,  0,-17,
        /* S */   0, -3,  0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2,  6,  0, -5, -7, -6, -4, -5,  0,-17,
        /* T */  -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -5, -3, -4, -9, -4,  0,  7,-13, -6, -3, -5, -4,  0,-17,
        /* W */ -13, -2, -8,-15,-15,-13,-17,-15, -7,-14,-10,-12,-13, -4,-14, -5,-13, 13, -5,-15, -7,-13,  0,-17,
        /* Y */  -8,-10, -4,-11, -4,-12, -8,-14, -3, -6, -7, -9,-11,  2,-13, -7, -6, -5, 10, -7,-10,-11,  0,-17,
        /* V */  -2, -8, -8, -8, -6, -7, -6, -5, -6,  2,  0, -9, -1, -8, -6, -6, -3,-15, -7,  7, -9, -8,  0,-17,
        /* B */  -7,  5, -4, -7,-11, -3, -7, -8, -4, -6, -7,  5, -3,-12, -5, -4, -5, -7,-10, -9,  5,  1,  0,-17,
        /* Z */  -6, -1, -2, -3,-14,  4, -2, -7, -3, -7, -7,  4, -3,-14, -5, -5, -4,-13,-11, -8,  1,  4,  0,-17,
        /* X */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,-17,
        /* * */ -17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,  1
      };

      return _data;
    }
  };

  template <>
  struct ScoringMatrixData_<int, AminoAcid, AdaptedIdentity>
  {
    enum
    {
      VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
      TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };
    static inline const int* getData()
    {
      static const int _data[TAB_SIZE] =
      {
        //      A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        /* A */ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* R */ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* N */ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* D */ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* C */ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* Q */ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* E */ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* G */ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* H */ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* I */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* L */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* K */ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* M */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* F */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* P */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* S */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -17,
        /* T */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -17,
        /* W */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -17,
        /* Y */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -17,
        /* V */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -17,
        /* B */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -17,
        /* Z */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -17,
        /* X */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17,
        /* * */ -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, 1
      };

      return _data;
    }
  };

} // namespace seqan


namespace OpenMS
{
  /**
    @brief Calculates a consensus from multiple ID runs based on PEPs and sequence similarities.

    @note The similarity scoring is based on an amino acid substitution matrix. Therefore only the raw amino acid sequences, without post-translational modifications (PTMs), can be considered for similarity scoring - PTMs are ignored during this step. However, PTMs on peptides are retained and separate results are produced for differently-modified peptides.

    @htmlinclude OpenMS_ConsensusIDAlgorithmPEPMatrix.parameters
    
    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusIDAlgorithmPEPMatrix :
    public ConsensusIDAlgorithmSimilarity
  {
  public:
    /// Default constructor
    ConsensusIDAlgorithmPEPMatrix();

  private:
    /// SeqAn similarity scoring
    typedef ::seqan::Score<int, ::seqan::ScoreMatrix< ::seqan::AminoAcid, ::seqan::Default> > SeqAnScore;

    /// SeqAn amino acid sequence
    typedef ::seqan::String< ::seqan::AminoAcid> SeqAnSequence;

    /// Similarity scoring method
    SeqAnScore scoring_method_;

    /// Alignment data structure
    ::seqan::Align<SeqAnSequence, ::seqan::ArrayGaps> alignment_;

    /// Not implemented
    ConsensusIDAlgorithmPEPMatrix(const ConsensusIDAlgorithmPEPMatrix&);

    /// Not implemented
    ConsensusIDAlgorithmPEPMatrix& operator=(const ConsensusIDAlgorithmPEPMatrix&);

    /// Docu in base class
    void updateMembers_() override;

    /// Sequence similarity based on substitution matrix (ignores PTMs)
    double getSimilarity_(AASequence seq1, AASequence seq2) override;

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSIDALGORITHMPEPMATRIX_H
