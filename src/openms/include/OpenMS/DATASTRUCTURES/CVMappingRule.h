// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  class CVMappingTerm;

  /**
      @brief Representation of a CV Mapping rule used by CVMappings

      Representation of a controlled vocabulary mapping rule.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI CVMappingRule
  {
public:

    /// enum to specify the requirement level
    enum RequirementLevel
    {
      MUST = 0,
      SHOULD = 1,
      MAY = 2
    };

    /// enum to specify the combination operator
    enum CombinationsLogic
    {
      OR = 0,
      AND = 1,
      XOR = 2
    };

    /// Default constructor
    CVMappingRule();

    /// Copy constructor
    CVMappingRule(const CVMappingRule& rhs);

    /// Destructor
    virtual ~CVMappingRule();

    /// Assignment operator
    CVMappingRule& operator=(const CVMappingRule& rhs);

    /** @name Accessors
    */
    //@{
    /// sets the identifier of the rule
    void setIdentifier(const String& identifier);

    /// returns the identifier of the rule
    const String& getIdentifier() const;

    /// sets the path of the element, where this rule is allowed
    void setElementPath(const String& element_path);

    /// returns the path of the element, where this rule is allowed
    const String& getElementPath() const;

    /// sets the requirement level of this rule
    void setRequirementLevel(RequirementLevel level);

    /// returns the requirement level of this rule
    RequirementLevel getRequirementLevel() const;

    /// sets the combination operator of the rule
    void setCombinationsLogic(CombinationsLogic combinations_logic);

    /// returns the combinations operator of the rule
    CombinationsLogic getCombinationsLogic() const;

    /// sets the scope path of the rule
    void setScopePath(const String& path);

    /// returns the scope path of the rule
    const String& getScopePath() const;

    /// sets the terms which are allowed
    void setCVTerms(const std::vector<CVMappingTerm>& cv_terms);

    /// returns the allowed terms
    const std::vector<CVMappingTerm>& getCVTerms() const;

    /// adds a term to the allowed terms
    void addCVTerm(const CVMappingTerm& cv_terms);
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVMappingRule& rhs) const;

    /// inequality operator
    bool operator!=(const CVMappingRule& rhs) const;
    //@}

protected:

    String identifier_;

    String element_path_;

    RequirementLevel requirement_level_;

    String scope_path_;

    CombinationsLogic combinations_logic_;

    std::vector<CVMappingTerm> cv_terms_;
  };

} // namespace OpenMS

