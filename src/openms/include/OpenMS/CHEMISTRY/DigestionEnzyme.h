// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <functional>
#include <iosfwd>
#include <set>
#include <vector>

namespace OpenMS
{
  /**
     @ingroup Chemistry

     @brief Base class for digestion enzymes
  */
  class OPENMS_DLLAPI DigestionEnzyme
  {
  public:

    /** @name Constructors
    */
    //@{
    /// Copy constructor
    DigestionEnzyme(const DigestionEnzyme&) = default;

    /// Move constructor
    DigestionEnzyme(DigestionEnzyme&&) = default;

    /// Detailed constructor
    explicit DigestionEnzyme(const std::string& name,
                             const std::string& cleavage_regex,
                             const std::set<std::string>& synonyms = std::set<std::string>(),
                             std::string regex_description = "");

    /// Detailed constructor 2
    explicit DigestionEnzyme(const std::string& name,
                             std::string cut_before,
                             const std::string& nocut_after = "",
                             std::string sense = "C",
                             const std::set<std::string>& synonyms = std::set<std::string>(),
                             std::string regex_description = "");

    /// Destructor
    virtual ~DigestionEnzyme();
    //@}

    /** @name Assignment
     */
    //@{
    /// Assignment operator
    DigestionEnzyme& operator=(const DigestionEnzyme&) = default;

    /// Move assignment operator
    DigestionEnzyme& operator=(DigestionEnzyme&&) & = default;
    //@}

    /** Accessors
    */
    //@{
    /// sets the name of the enzyme
    void setName(const std::string& name);

    /// returns the name of the enzyme
    const std::string& getName() const;

    /// sets the synonyms
    void setSynonyms(const std::set<std::string>& synonyms);

    /// adds a synonym
    void addSynonym(const std::string& synonym);

    /// returns the synonyms
    const std::set<std::string>& getSynonyms() const;

    /// sets the cleavage regex
    void setRegEx(const std::string& cleavage_regex);

    /// returns the cleavage regex
    const std::string& getRegEx() const;

    /// sets the regex description
    void setRegExDescription(const std::string& value);

    /// returns the regex description
    const std::string& getRegExDescription() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const DigestionEnzyme& enzyme) const;

    /// inequality operator
    bool operator!=(const DigestionEnzyme& enzyme) const;

    /// equality operator for regex
    bool operator==(const std::string& cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(const std::string& cleavage_regex) const;

    /// order operator
    bool operator<(const DigestionEnzyme& enzyme) const;
    //@}

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    virtual bool setValueFromFile(const std::string& key, const std::string& value);

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme);

  protected:

    /// default constructor
    DigestionEnzyme();

    // basic
    std::string name_;

    std::string cleavage_regex_;

    std::set<std::string> synonyms_;

    std::string regex_description_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme);
} // namespace OpenMS

namespace std
{
  /// std::hash specialization for DigestionEnzyme
  template<>
  struct hash<OpenMS::DigestionEnzyme>
  {
    std::size_t operator()(const OpenMS::DigestionEnzyme& enzyme) const noexcept
    {
      std::size_t seed = OpenMS::fnv1a_hash_string(enzyme.getName());
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(enzyme.getRegEx()));
      OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(enzyme.getRegExDescription()));
      // Hash the synonyms set
      for (const auto& syn : enzyme.getSynonyms())
      {
        OpenMS::hash_combine(seed, OpenMS::fnv1a_hash_string(syn));
      }
      return seed;
    }
  };
} // namespace std

