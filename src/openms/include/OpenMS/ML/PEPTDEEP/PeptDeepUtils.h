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

namespace OpenMS {
  namespace ML {

    // --- Shared PeptDeep Architecture Constants ---
    constexpr int64_t PEPTDEEP_MOD_ELEMENTS = 109;

    /// @brief AlphaPeptDeep's exact 109-element array for modification tensor mapping.
    /// Elements missing from this list (e.g., isotopes like "13C") are gracefully binned into the final "Other" channel.
    const std::vector<std::string> ALPHAPEPTDEEP_MOD_ELEMENTS = {
        "C", "H", "N", "O", "P", "S", "Br", "Cl", "F", "Fe", "I", "K", "Na", "Zn",
        "Se", "Mg", "Ca", "Cu", "Mn", "Ni", "Mo", "Ag", "Co", "Au", "V", "Pt",
        "Ru", "Cd", "Cr", "W", "Pb", "Li", "Rb", "Cs", "Fr", "Be", "Sr", "Ba",
        "Ra", "Sc", "Y", "Ti", "Zr", "Hf", "Rf", "Nb", "Ta", "Db", "Tc", "Re",
        "Bh", "Os", "Hs", "Rh", "Ir", "Mt", "Pd", "Ds", "Rg", "Cn", "Nh", "Fl",
        "Mc", "Lv", "Ts", "Og", "B", "Al", "Ga", "In", "Tl", "Si", "Ge", "Sn",
        "Uut", "As", "Sb", "Bi", "Uup", "Te", "Po", "At", "Uus", "Rn", "Uuo",
        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
        "Yb", "Lu", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
        "Fm", "Md", "No", "Lr"
    };

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