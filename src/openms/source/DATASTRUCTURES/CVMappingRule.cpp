// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>

using namespace std;

namespace OpenMS
{
  // CV mapping rule implementation
  CVMappingRule::CVMappingRule() :
    requirement_level_(CVMappingRule::MUST),
    combinations_logic_(CVMappingRule::OR)
  {
  }

  CVMappingRule::CVMappingRule(const CVMappingRule& rhs) = default;

  CVMappingRule::~CVMappingRule() = default;

  CVMappingRule& CVMappingRule::operator=(const CVMappingRule& rhs)
  {
    if (this != &rhs)
    {
      identifier_ = rhs.identifier_;
      element_path_ = rhs.element_path_;
      requirement_level_ = rhs.requirement_level_;
      scope_path_ = rhs.scope_path_;
      combinations_logic_ = rhs.combinations_logic_;
      cv_terms_ = rhs.cv_terms_;
    }
    return *this;
  }

  bool CVMappingRule::operator==(const CVMappingRule& rhs) const
  {
    return identifier_ == rhs.identifier_ &&
           element_path_ == rhs.element_path_ &&
           requirement_level_ == rhs.requirement_level_ &&
           scope_path_ == rhs.scope_path_ &&
           combinations_logic_ == rhs.combinations_logic_ &&
           cv_terms_ == rhs.cv_terms_;
  }

  bool CVMappingRule::operator!=(const CVMappingRule& rhs) const
  {
    return !(*this == rhs);
  }

  void CVMappingRule::setIdentifier(const String& identifier)
  {
    identifier_ = identifier;
  }

  const String& CVMappingRule::getIdentifier() const
  {
    return identifier_;
  }

  void CVMappingRule::setElementPath(const String& element_path)
  {
    element_path_ = element_path;
  }

  const String& CVMappingRule::getElementPath() const
  {
    return element_path_;
  }

  void CVMappingRule::setRequirementLevel(RequirementLevel level)
  {
    requirement_level_ = level;
  }

  CVMappingRule::RequirementLevel CVMappingRule::getRequirementLevel() const
  {
    return requirement_level_;
  }

  void CVMappingRule::setCombinationsLogic(CombinationsLogic combinations_logic)
  {
    combinations_logic_ = combinations_logic;
  }

  CVMappingRule::CombinationsLogic CVMappingRule::getCombinationsLogic() const
  {
    return combinations_logic_;
  }

  void CVMappingRule::setScopePath(const String& path)
  {
    scope_path_ = path;
  }

  const String& CVMappingRule::getScopePath() const
  {
    return scope_path_;
  }

  void CVMappingRule::setCVTerms(const vector<CVMappingTerm>& cv_terms)
  {
    cv_terms_ = cv_terms;
  }

  const vector<CVMappingTerm>& CVMappingRule::getCVTerms() const
  {
    return cv_terms_;
  }

  void CVMappingRule::addCVTerm(const CVMappingTerm& cv_term)
  {
    cv_terms_.push_back(cv_term);
  }

} // namespace OpenMS
