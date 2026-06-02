// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/CVTermListInterface.h>
#include <OpenMS/CONCEPT/HashUtils.h>

#include <cstddef>

namespace OpenMS
{

  /// Helper template to hash any type with getCVTerms() method (CVTermList, CVTermListInterface)
  template<typename T>
  inline std::size_t hashCVTerms(const T& obj) noexcept
  {
    std::size_t seed = 0;
    const auto& cv_terms = obj.getCVTerms();
    for (const auto& [accession, terms] : cv_terms)
    {
      hash_combine(seed, fnv1a_hash_string(accession));
      for (const auto& term : terms)
      {
        hash_combine(seed, fnv1a_hash_string(term.getAccession()));
        hash_combine(seed, fnv1a_hash_string(term.getName()));
        hash_combine(seed, fnv1a_hash_string(term.getCVIdentifierRef()));
        if (term.hasValue())
        {
          hash_combine(seed, fnv1a_hash_string(term.getValue().toString()));
        }
        if (term.hasUnit())
        {
          hash_combine(seed, fnv1a_hash_string(term.getUnit().accession));
        }
      }
    }
    return seed;
  }

  // Convenience wrappers for backward compatibility
  inline std::size_t hashCVTermList(const CVTermList& cvtl) noexcept { return hashCVTerms(cvtl); }
  inline std::size_t hashCVTermListInterface(const CVTermListInterface& cvtli) noexcept { return hashCVTerms(cvtli); }

} // namespace OpenMS
