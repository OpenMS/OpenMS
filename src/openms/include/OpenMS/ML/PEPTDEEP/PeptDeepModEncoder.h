// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <cstddef>
#include <string>
#include <vector>

namespace OpenMS
{
  class AASequence;

  namespace ML
  {
    /// @brief Encodes OpenMS modifications into PeptDeep's mod_x tensor.
    ///
    /// PeptDeep describes a modification purely by the elemental composition of its
    /// mass delta: the per-site feature vector holds the atom counts of the
    /// modification's diff formula, laid out over a fixed 109-entry element list
    /// (AlphaPeptDeep's @c mod_elements, ending in the isotopes 2H/13C/15N/18O and a
    /// catch-all "?" slot). Nothing about the modification's name or UniMod identity
    /// enters the model, so any OpenMS modification with a diff formula can be encoded
    /// without a name-to-name mapping table.
    ///
    /// Row layout matches AlphaPeptDeep's @c parse_mod_feature: a peptide of length
    /// @c nAA occupies @c nAA+2 rows, where row 0 is the N-terminus, rows 1..nAA are
    /// the residues, and row @c nAA+1 is the C-terminus.
    class OPENMS_DLLAPI PeptDeepModEncoder
    {
    public:
      /// Length of the per-site modification feature vector (AlphaPeptDeep mod_feature_size).
      static constexpr size_t MOD_FEATURE_SIZE = 109;

      /// @brief Index of an element symbol in PeptDeep's mod feature vector.
      /// @param symbol Element symbol as used by OpenMS, e.g. "C" or the isotope form "(13)C".
      /// @return The element's index, or the index of the catch-all "?" entry for
      ///         anything PeptDeep does not model explicitly.
      static size_t elementIndex(const std::string& symbol);

      /// @brief Index of the catch-all "?" entry.
      static size_t unknownElementIndex();

      /// @brief Adds one peptide's modifications into a mod_x block.
      /// @param seq The (possibly modified) peptide.
      /// @param mod_x Pointer to the first of @c rows * MOD_FEATURE_SIZE floats for this
      ///        peptide, laid out row-major. The caller owns the buffer and must have
      ///        zero-initialised it; this function only adds.
      /// @param rows Number of rows available, i.e. the padded sequence length. Must be at
      ///        least @c seq.size()+2; surplus rows (padding) are left untouched.
      /// @throws Exception::InvalidValue if @c rows is too small for the peptide.
      ///
      /// Contributions are accumulated, so several modifications on one site add up, as
      /// they do in AlphaPeptDeep.
      static void encode(const AASequence& seq, float* mod_x, size_t rows);

      /// @brief The element symbols behind the feature vector, in model order.
      /// Exposed for testing and diagnostics.
      static const std::vector<std::string>& elementOrder();
    };
  } // namespace ML
} // namespace OpenMS
