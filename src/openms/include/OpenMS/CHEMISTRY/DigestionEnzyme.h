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
#include <OpenMS/DATASTRUCTURES/String.h>

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
    explicit DigestionEnzyme(const String& name,
                             const String& cleavage_regex,
                             const std::set<String>& synonyms = std::set<String>(),
                             String regex_description = "");

    /// Detailed constructor 2
    explicit DigestionEnzyme(const String& name,
                             String cut_before,
                             const String& nocut_after = "",
                             String sense = "C",
                             const std::set<String>& synonyms = std::set<String>(),
                             String regex_description = "");

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
    void setName(const String& name);

    /// returns the name of the enzyme
    const String& getName() const;

    /// sets the synonyms
    void setSynonyms(const std::set<String>& synonyms);

    /// adds a synonym
    void addSynonym(const String& synonym);

    /// returns the synonyms
    const std::set<String>& getSynonyms() const;

    /// sets the cleavage regex
    void setRegEx(const String& cleavage_regex);

    /// returns the cleavage regex
    const String& getRegEx() const;

    /// sets the regex description
    void setRegExDescription(const String& value);

    /// returns the regex description
    const String& getRegExDescription() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const DigestionEnzyme& enzyme) const;

    /// inequality operator
    bool operator!=(const DigestionEnzyme& enzyme) const;

    /// equality operator for regex
    bool operator==(const String& cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(const String& cleavage_regex) const;

    /// order operator
    bool operator<(const DigestionEnzyme& enzyme) const;
    //@}

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    virtual bool setValueFromFile(const String& key, const String& value);

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme);

  protected:

    /// default constructor
    DigestionEnzyme();

    // basic
    String name_;

    String cleavage_regex_;

    std::set<String> synonyms_;

    String regex_description_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzyme& enzyme);
}

