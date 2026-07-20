// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
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

    /**
     * @brief Generates an empty mod_x tensor for unmodified peptides.
     * Resolves to a flat vector of zeros of size (batch_size * sequence_length * 109).
     */
    inline std::vector<float> generateUnmodifiedModXTensor(size_t batch_size, size_t sequence_length) {
        return std::vector<float>(batch_size * sequence_length * PEPTDEEP_MOD_ELEMENTS, 0.0f);
    }

  }
}