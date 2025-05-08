// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{
  ModificationDefinition::ModificationDefinition() :
    mod_(nullptr),
    fixed_modification_(true),
    max_occurrences_(0)
  {
  }

  ModificationDefinition::ModificationDefinition(const ModificationDefinition& rhs) = default;

  ModificationDefinition::ModificationDefinition(const String& mod, bool fixed, UInt max_occur) :
    mod_(nullptr),
    fixed_modification_(fixed),
    max_occurrences_(max_occur)
  {
    setModification(mod);
  }

  ModificationDefinition::ModificationDefinition(const ResidueModification& mod, bool fixed, UInt max_occur) :
    mod_(&mod),
    fixed_modification_(fixed),
    max_occurrences_(max_occur)
  {
  }

  ModificationDefinition& ModificationDefinition::operator=(const ModificationDefinition& rhs)
  {
    if (this != &rhs)
    {
      mod_ = rhs.mod_;
      fixed_modification_ = rhs.fixed_modification_;
      max_occurrences_ = rhs.max_occurrences_;
    }
    return *this;
  }

  bool ModificationDefinition::operator==(const ModificationDefinition& rhs) const
  {
    return mod_ == rhs.mod_ &&
           fixed_modification_ == rhs.fixed_modification_ &&
           max_occurrences_ == rhs.max_occurrences_;
  }

  bool ModificationDefinition::operator!=(const ModificationDefinition& rhs) const
  {
    return !(*this == rhs);
  }

  ModificationDefinition::~ModificationDefinition() = default;

  bool ModificationDefinition::operator<(const ModificationDefinition& rhs) const
  {
    return this->getModificationName() < rhs.getModificationName();
  }

  void ModificationDefinition::setFixedModification(bool fixed_mod)
  {
    fixed_modification_ = fixed_mod;
  }

  bool ModificationDefinition::isFixedModification() const
  {
    return fixed_modification_;
  }

  void ModificationDefinition::setModification(const String& modification)
  {
    mod_ = ModificationsDB::getInstance()->getModification(modification);
  }

  const ResidueModification& ModificationDefinition::getModification() const
  {
    if (!mod_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No modification defined", nullptr);
    }
    return *mod_;
  }

  String ModificationDefinition::getModificationName() const
  {
    if (mod_ != nullptr)
    {
      return mod_->getFullId();
    }
    return "";
  }

  void ModificationDefinition::setMaxOccurrences(UInt max_occurrences)
  {
    max_occurrences_ = max_occurrences;
  }

  UInt ModificationDefinition::getMaxOccurrences() const
  {
    return max_occurrences_;
  }

} // namespace OpenMS
