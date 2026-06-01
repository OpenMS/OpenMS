// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <cstdint>
#include <functional>

namespace OpenMS
{
  /** @ingroup Chemistry

          @brief Representation of modification definition

          This class defines a modification type e.g. a input parameter of a search engine.
          The modification is defined using an unique name of the modification present
          in the modifications DB instance.
  */
  class OPENMS_DLLAPI ModificationDefinition
  {
public:

    /** @name Constructor and Destructors
    */
    //@{
    /// default constructor
    ModificationDefinition();

    /// copy constructor
    ModificationDefinition(const ModificationDefinition& rhs);

    /// detailed constructor specifying the modification by name
    explicit ModificationDefinition(const String& mod, bool fixed = true, UInt max_occur = 0);

    /// direct constructor from a residue modification
    explicit ModificationDefinition(const ResidueModification& mod, bool fixed = true, UInt max_occur = 0);

    /// destructor
    virtual ~ModificationDefinition();
    //@}

    /** @name Accessors
    */
    //@{
    /// sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
    void setFixedModification(bool fixed);

    /// returns if the modification if fixed true, else false
    bool isFixedModification() const;

    /// set the maximal number of occurrences per peptide (unbounded if 0)
    void setMaxOccurrences(UInt num);

    /// returns the maximal number of occurrences per peptide
    UInt getMaxOccurrences() const;

    /// returns the name of the modification
    String getModificationName() const;

    /// sets the modification, allowed are unique names provided by ModificationsDB
    void setModification(const String& modification);

    /**
       @brief Returns the modification

       @throw Exception::InvalidValue if no modification was set
    */
    const ResidueModification& getModification() const;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    ModificationDefinition& operator=(const ModificationDefinition& element);
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const ModificationDefinition& rhs) const;

    /// inequality operator
    bool operator!=(const ModificationDefinition& rhs) const;

    /// less than operator for e.g. usage in maps; only mod FullIds are compared!
    bool operator<(const OpenMS::ModificationDefinition&) const;
    //@}

protected:

    /// the modification
    const ResidueModification* mod_;

    /// fixed (true) or variable (false)
    bool fixed_modification_;

    /// maximal number of occurrences per peptide
    UInt max_occurrences_;
  };

} // namespace OpenMS

namespace std
{
  /**
   * @brief Hash function for OpenMS::ModificationDefinition.
   *
   * Computes a hash based on all fields used in operator==:
   * - mod_ pointer (pointer identity, not content)
   * - fixed_modification_ (bool)
   * - max_occurrences_ (UInt)
   *
   * @note Since operator== compares pointer addresses (not dereferencing),
   *       this hash uses pointer identity. Hash is consistent within a process
   *       but not reproducible across process runs due to address space layout.
   */
  template<>
  struct hash<OpenMS::ModificationDefinition>
  {
    std::size_t operator()(const OpenMS::ModificationDefinition& md) const noexcept
    {
      // Hash the modification pointer as uintptr_t (matches pointer comparison in operator==)
      std::size_t seed = OpenMS::hash_int(reinterpret_cast<std::uintptr_t>(&md.getModification()));
      // Hash fixed_modification_ flag
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(md.isFixedModification())));
      // Hash max_occurrences_
      OpenMS::hash_combine(seed, OpenMS::hash_int(md.getMaxOccurrences()));
      return seed;
    }
  };
} // namespace std

