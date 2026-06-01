// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <functional>

namespace OpenMS
{
  /**
      @brief Description of the software used for processing

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Software :
    public CVTermList
  {
public:
    /// Constructor
    explicit Software(const String& name = "", const String& version = "");
    /// Copy constructor
    Software(const Software&) = default;
    /// Move constructor
    Software(Software&&) = default;
    /// Destructor
    ~Software() override;

    /// Assignment operator
    Software& operator=(const Software&) = default;
    /// Move assignment operator
    Software& operator=(Software&&)& = default;

    /// Equality operator
    bool operator==(const Software& rhs) const;
    /// Inequality operator
    bool operator!=(const Software& rhs) const;
    /// Less-than operator (for sorting)
    bool operator<(const Software& rhs) const;

    /// Returns the name of the software
    const String& getName() const;
    /// Sets the name of the software
    void setName(const String& name);

    /// Returns the software version
    const String& getVersion() const;
    /// Sets the software version
    void setVersion(const String& version);

protected:
    String name_;
    String version_;
  };
} // namespace OpenMS

// Hash function specialization for Software
namespace std
{
  /**
   * @brief Hash function for OpenMS::Software.
   *
   * Hashes name and version fields.
   * Note: CVTermList base class fields are not included in the hash as CVTermList
   * does not yet have a hash implementation. This hash is consistent with operator==
   * for objects that differ only in the Software-specific fields.
   */
  template<>
  struct hash<OpenMS::Software>
  {
    std::size_t operator()(const OpenMS::Software& s) const noexcept
    {
      std::size_t seed = OpenMS::fnv1a_hash_string(s.getName());
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(s.getVersion()));
      return seed;
    }
  };
} // namespace std

