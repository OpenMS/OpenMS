// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>

using namespace std;

namespace OpenMS
{
  // CV term implementation
  CVMappingTerm::CVMappingTerm() :
    use_term_name_(false),
    use_term_(false),
    is_repeatable_(false),
    allow_children_(false)
  {
  }

  CVMappingTerm::CVMappingTerm(const CVMappingTerm& rhs) = default;

  CVMappingTerm::~CVMappingTerm() = default;

  CVMappingTerm& CVMappingTerm::operator=(const CVMappingTerm& rhs)
  {
    if (this != &rhs)
    {
      accession_ = rhs.accession_;
      use_term_name_ = rhs.use_term_name_;
      use_term_ = rhs.use_term_;
      term_name_ = rhs.term_name_;
      is_repeatable_ = rhs.is_repeatable_;
      allow_children_ = rhs.allow_children_;
      cv_identifier_ref_ = rhs.cv_identifier_ref_;
    }
    return *this;
  }

  bool CVMappingTerm::operator==(const CVMappingTerm& rhs) const
  {
    return accession_ == rhs.accession_ &&
           use_term_name_ == rhs.use_term_name_ &&
           use_term_ == rhs.use_term_ &&
           term_name_ == rhs.term_name_ &&
           is_repeatable_ == rhs.is_repeatable_ &&
           allow_children_ == rhs.allow_children_ &&
           cv_identifier_ref_ == rhs.cv_identifier_ref_;
  }

  bool CVMappingTerm::operator!=(const CVMappingTerm& rhs) const
  {
    return !(*this == rhs);
  }

  void CVMappingTerm::setAccession(const String& accession)
  {
    accession_ = accession;
  }

  const String& CVMappingTerm::getAccession() const
  {
    return accession_;
  }

  void CVMappingTerm::setUseTermName(bool use_term_name)
  {
    use_term_name_ = use_term_name;
  }

  bool CVMappingTerm::getUseTermName() const
  {
    return use_term_name_;
  }

  void CVMappingTerm::setUseTerm(bool use_term)
  {
    use_term_ = use_term;
  }

  bool CVMappingTerm::getUseTerm() const
  {
    return use_term_;
  }

  void CVMappingTerm::setTermName(const String& term_name)
  {
    term_name_ = term_name;
  }

  const String& CVMappingTerm::getTermName() const
  {
    return term_name_;
  }

  void CVMappingTerm::setIsRepeatable(bool is_repeatable)
  {
    is_repeatable_ = is_repeatable;
  }

  bool CVMappingTerm::getIsRepeatable() const
  {
    return is_repeatable_;
  }

  void CVMappingTerm::setAllowChildren(bool allow_children)
  {
    allow_children_ = allow_children;
  }

  bool CVMappingTerm::getAllowChildren() const
  {
    return allow_children_;
  }

  void CVMappingTerm::setCVIdentifierRef(const String& cv_identifier_ref)
  {
    cv_identifier_ref_ = cv_identifier_ref;
  }

  const String& CVMappingTerm::getCVIdentifierRef() const
  {
    return cv_identifier_ref_;
  }

} // namespace OpenMS
