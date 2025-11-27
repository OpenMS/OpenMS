// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTermList.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace std;

namespace OpenMS
{

  CVTermList::~CVTermList() = default;

  CVTermList::CVTermList(CVTermList&& rhs) noexcept :
    MetaInfoInterface(std::move(rhs)),
    cv_terms_(std::move(rhs.cv_terms_))
  {
  }

  void CVTermList::addCVTerm(const CVTerm& cv_term)
  {
    // TODO exception if empty
    cv_terms_[cv_term.getAccession()].push_back(cv_term);
  }

  void CVTermList::setCVTerms(const vector<CVTerm>& cv_terms)
  {
    for (const CVTerm& tr : cv_terms)
    {
      addCVTerm(tr);
    }
    return;
  }

  void CVTermList::replaceCVTerm(const CVTerm& cv_term)
  {
    std::vector<CVTerm> tmp;
    tmp.push_back(cv_term);
    cv_terms_[cv_term.getAccession()] = tmp;
  }

  void CVTermList::replaceCVTerms(const vector<CVTerm>& cv_terms, const String& accession)
  {
    cv_terms_[accession] = cv_terms;
  }

  void CVTermList::replaceCVTerms(const std::map<String, vector<CVTerm> >& cv_term_map)
  {
    cv_terms_ = cv_term_map;
  }

  void CVTermList::consumeCVTerms(const std::map<String, vector<CVTerm> >& cv_term_map)
  {
    for (std::map<String, std::vector<CVTerm> >::const_iterator it = cv_term_map.begin(); it != cv_term_map.end(); ++it)
    {
      cv_terms_[it->first].insert(cv_terms_[it->first].end(), it->second.begin(), it->second.end());
    }
  }

  const std::map<String, vector<CVTerm> >& CVTermList::getCVTerms() const
  {
    return cv_terms_;
  }

  bool CVTermList::hasCVTerm(const String& accession) const
  {
    return cv_terms_.find(accession) != cv_terms_.end();
  }

  bool CVTermList::operator==(const CVTermList& cv_term_list) const
  {
    return MetaInfoInterface::operator==(cv_term_list) && cv_terms_ == cv_term_list.cv_terms_;
  }

  bool CVTermList::operator!=(const CVTermList& cv_term_list) const
  {
    return !(*this == cv_term_list);
  }

  bool CVTermList::empty() const
  {
    return cv_terms_.empty();
  }

} // namespace OpenMS

