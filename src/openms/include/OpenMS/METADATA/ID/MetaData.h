// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <cstdint> // for "uintptr_t"

namespace OpenMS
{
  namespace IdentificationDataInternal
  {
    /// Wrapper that adds @p operator< to iterators, so they can be used as (part of) keys in maps/sets or @p multi_index_containers
    template <typename Iterator>
    struct IteratorWrapper: public Iterator
    {
      IteratorWrapper(): Iterator() {}

      IteratorWrapper(const Iterator& it): Iterator(it) {}

      bool operator<(const IteratorWrapper& other) const
      {
        // compare by address of referenced element:
        return &(**this) < &(*other);
      }

      /// Conversion to pointer type for hashing
      operator uintptr_t() const
      {
        return uintptr_t(&(**this));
      }
    };


    enum MoleculeType
    {
      PROTEIN,
      COMPOUND,
      RNA
    };


    enum MassType
    {
      MONOISOTOPIC,
      AVERAGE
    };
  }
}
