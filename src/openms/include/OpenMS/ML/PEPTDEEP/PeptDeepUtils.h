// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav, Justin Sing $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <string>
#include <vector>
#include <array>
#include <string_view>

namespace OpenMS {
  namespace ML {

    // --- Shared PeptDeep Architecture Constants ---
    constexpr int64_t PEPTDEEP_MOD_ELEMENTS = 109;

    /// @brief AlphaPeptDeep's exact 109-element array for modification tensor mapping.
    /// Elements missing from this list are binned into the final "Other" channel ("?").
    /// Heavy isotopes (2H, 13C, 15N, 18O) map strictly to their dedicated channels.
    inline constexpr std::array<std::string_view, 109> ALPHAPEPTDEEP_MOD_ELEMENTS = {
    "C", "H", "N", "O", "P", "S", "B", "F", "I", "K", "U", "V", "W", "X", "Y",
    "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "Ba", "Be", "Bi", "Bk",
    "Br", "Ca", "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Dy",
    "Er", "Es", "Eu", "Fe", "Fm", "Fr", "Ga", "Gd", "Ge", "He", "Hf", "Hg",
    "Ho", "In", "Ir", "Kr", "La", "Li", "Lr", "Lu", "Md", "Mg", "Mn", "Mo",
    "Na", "Nb", "Nd", "Ne", "Ni", "No", "Np", "Os", "Pa", "Pb", "Pd", "Pm",
    "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rh", "Rn", "Ru", "Sb", "Sc",
    "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Th", "Ti", "Tl",
    "Tm", "Xe", "Yb", "Zn", "Zr", "2H", "13C", "15N", "18O", "?"
    };

    /**
     * @brief Width of the instrument one-hot in AlphaPeptDeep's MS2 model.
     *
     * Mirrors @p max_instrument_num in peptdeep's model_const.yaml. Valid instrument indices are
     * [0, PEPTDEEP_MAX_INSTRUMENT_NUM): the MS2 graph feeds the index through a OneHot of this
     * width, so the bound is a property of the trained weights rather than a choice. Out of range
     * is not an error to ONNX -- its OneHot answers with an all-off row, which silently drops the
     * instrument from the prediction -- so it has to be rejected here instead.
     * tools/scripts/export_peptdeep_models_to_onnx.py checks an exported graph against this.
     */
    constexpr int64_t PEPTDEEP_MAX_INSTRUMENT_NUM = 8;

    /// @brief Index peptdeep gives an instrument it does not recognise (featurize.py: max_instrument_num - 1).
    constexpr int64_t PEPTDEEP_UNKNOWN_INSTRUMENT_INDEX = PEPTDEEP_MAX_INSTRUMENT_NUM - 1;

    /**
     * @brief AlphaPeptDeep's instrument names, in the order that *is* their encoding.
     *
     * From model_const.yaml, whose own comment reads "We MUST keep the order of these instruments
     * for models": position is the index the models were trained against, so entries may only be
     * appended, never reordered or removed. The slots between the last name and
     * PEPTDEEP_UNKNOWN_INSTRUMENT_INDEX are the room peptdeep has left for new instruments.
     */
    inline constexpr std::array<std::string_view, 5> ALPHAPEPTDEEP_INSTRUMENTS = {
      "QE", "Lumos", "timsTOF", "SciexTOF", "ThermoTOF"
    };

    /**
     * @brief Instrument index for @p name, mirroring peptdeep's parse_instrument_indices().
     *
     * Case-insensitive, as peptdeep's own lookup is, and a name that is not in
     * ALPHAPEPTDEEP_INSTRUMENTS lands on the catch-all slot rather than on a neighbouring
     * instrument. Prefer this over writing index literals: the order is peptdeep's, not ours.
     *
     * @param name Instrument name, e.g. "QE" or "timsTOF".
     * @return Index in [0, PEPTDEEP_MAX_INSTRUMENT_NUM).
     */
    inline int64_t instrumentIndex(std::string_view name)
    {
      auto upper = [](char c) { return (c >= 'a' && c <= 'z') ? static_cast<char>(c - ('a' - 'A')) : c; };

      for (Size i = 0; i < ALPHAPEPTDEEP_INSTRUMENTS.size(); ++i)
      {
        const std::string_view candidate = ALPHAPEPTDEEP_INSTRUMENTS[i];
        if (candidate.size() != name.size()) { continue; }

        bool same = true;
        for (Size c = 0; c < name.size(); ++c)
        {
          if (upper(name[c]) != upper(candidate[c])) { same = false; break; }
        }
        if (same) { return static_cast<int64_t>(i); }
      }
      return PEPTDEEP_UNKNOWN_INSTRUMENT_INDEX;
    }

    /**
     * @brief Validates instrument indices for PeptDeep MS2 inference.
     *
     * @param instrument_indices One index per peptide.
     * @throws Exception::IllegalArgument if an index is outside the one-hot the model encodes.
     */
    inline void validateInstrumentIndices(const std::vector<int64_t>& instrument_indices)
    {
      for (int64_t index : instrument_indices)
      {
        if (index < 0 || index >= PEPTDEEP_MAX_INSTRUMENT_NUM)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Instrument index " + std::to_string(index) + " is outside the range [0, "
            + std::to_string(PEPTDEEP_MAX_INSTRUMENT_NUM) + ") that the MS2 model encodes. Derive "
            "indices from ML::instrumentIndex() rather than passing literals.");
        }
      }
    }

    /**
     * @brief Maps amino acid characters to 1-based token indices for PeptDeep models.
     *
     * Converts uppercase and lowercase letters A-Z to indices 1-26 using canonical
     * ord-offset encoding. Non-alphabetic characters map to 0, which serves as the
     * padding and unknown token in PeptDeep ONNX input tensors.
     *
     * @param aa Amino acid character (case-insensitive)
     * @return Token index: 1-26 for A-Z/a-z, 0 for padding/unknown
     */
    inline OpenMS::Int64 getAAIndex(char aa) {
      if (aa >= 'A' && aa <= 'Z') return aa - 'A' + 1;
      if (aa >= 'a' && aa <= 'z') return aa - 'a' + 1;
      return 0; // 0 serves as the padding and unknown token
    }

    /**
     * @brief Validates a peptide sequence for PeptDeep inference.
     * Complex character validation and modification parsing is handled downstream by AASequence.
     * @param peptide The raw peptide string to validate.
     */
    inline void validatePeptide(const std::string& peptide) {
        if (peptide.empty()) {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptide sequence cannot be empty.");
        }
    }

  }
}