// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS {
  namespace ML {
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
  }
}