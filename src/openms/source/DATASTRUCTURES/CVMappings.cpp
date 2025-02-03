// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>

#include <iostream>

namespace OpenMS
{
  class CVMappingRule;
} // namespace OpenMS

using namespace std;

namespace OpenMS
{
  CVMappings::CVMappings() = default;

  CVMappings::CVMappings(const CVMappings& rhs) = default;

  CVMappings::~CVMappings() = default;

  CVMappings& CVMappings::operator=(const CVMappings& rhs)
  {
    if (this != &rhs)
    {
      mapping_rules_ = rhs.mapping_rules_;
      cv_references_ = rhs.cv_references_;
      cv_references_vector_ = rhs.cv_references_vector_;
    }
    return *this;
  }

  bool CVMappings::operator==(const CVMappings& rhs) const
  {
    return mapping_rules_ == rhs.mapping_rules_ &&
           cv_references_ == rhs.cv_references_ &&
           cv_references_vector_ == rhs.cv_references_vector_;
  }

  bool CVMappings::operator!=(const CVMappings& rhs) const
  {
    return !(*this == rhs);
  }

  void CVMappings::setMappingRules(const vector<CVMappingRule>& cv_mapping_rules)
  {
    mapping_rules_ = cv_mapping_rules;
  }

  const vector<CVMappingRule>& CVMappings::getMappingRules() const
  {
    return mapping_rules_;
  }

  void CVMappings::addMappingRule(const CVMappingRule& cv_mapping_rule)
  {
    mapping_rules_.push_back(cv_mapping_rule);
  }

  void CVMappings::setCVReferences(const vector<CVReference>& cv_references)
  {
    for (vector<CVReference>::const_iterator it = cv_references.begin(); it != cv_references.end(); ++it)
    {
      cv_references_[it->getIdentifier()] = *it;
      cv_references_vector_.push_back(*it);
    }
  }

  const vector<CVReference>& CVMappings::getCVReferences() const
  {
    return cv_references_vector_;
  }

  void CVMappings::addCVReference(const CVReference& cv_reference)
  {
    if (hasCVReference(cv_reference.getIdentifier()))
    {
      cerr << "CVMappings: Warning: CV reference with identifier '" << cv_reference.getIdentifier() << "' already existing, ignoring it!" << endl;
      return;
    }
    cv_references_[cv_reference.getIdentifier()] = cv_reference;
    cv_references_vector_.push_back(cv_reference);
  }

  bool CVMappings::hasCVReference(const String& identifier)
  {
    return cv_references_.find(identifier) != cv_references_.end();
  }

} // namespace OpenMS
