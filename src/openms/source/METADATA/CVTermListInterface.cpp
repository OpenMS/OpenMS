// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTermListInterface.h>

#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/Helpers.h>
#include <map>

namespace OpenMS
{

  static const std::map<String, std::vector<CVTerm> > empty_cvterm_map = std::map<String, std::vector<CVTerm> >();

  CVTermListInterface::CVTermListInterface() :
      MetaInfoInterface(),
      cvt_ptr_(nullptr)
    {}

  CVTermListInterface::CVTermListInterface(const CVTermListInterface & rhs) :
    MetaInfoInterface(rhs),
    cvt_ptr_(nullptr)
  {
    if (rhs.cvt_ptr_ != nullptr)
    {
      cvt_ptr_ = new CVTermList(*rhs.cvt_ptr_);
    }
  }

  /// Move constructor
  CVTermListInterface::CVTermListInterface(CVTermListInterface&& rhs) noexcept :
    MetaInfoInterface(std::move(rhs)), // NOTE: rhs itself is an lvalue
    cvt_ptr_(rhs.cvt_ptr_)
  {
    // see http://thbecker.net/articles/rvalue_references/section_05.html
    // take ownership
    rhs.cvt_ptr_ = nullptr;
  }

  CVTermListInterface::~CVTermListInterface() 
  {
    delete cvt_ptr_;
  }

  CVTermListInterface & CVTermListInterface::operator=(const CVTermListInterface & rhs)
  {
    if (this != &rhs)
    {
      MetaInfoInterface::operator=(rhs);

      delete cvt_ptr_;
      cvt_ptr_ = nullptr;
      if (rhs.cvt_ptr_ != nullptr)
      {
        cvt_ptr_ = new CVTermList(*rhs.cvt_ptr_);
      }
    }
    return *this;
  }

  CVTermListInterface& CVTermListInterface::operator=(CVTermListInterface&& rhs) noexcept
  {
    if (&rhs == this)
    {
      return *this;
    }

    MetaInfoInterface::operator=(std::move(rhs));

    // free memory and assign rhs memory
    delete cvt_ptr_;
    cvt_ptr_ = rhs.cvt_ptr_;
    rhs.cvt_ptr_ = nullptr;

    return *this;
  }

  bool CVTermListInterface::operator==(const CVTermListInterface& rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           Helpers::cmpPtrSafe<CVTermList*>(cvt_ptr_, rhs.cvt_ptr_);
  }

  bool CVTermListInterface::operator!=(const CVTermListInterface& rhs) const
  {
    return !(*this == rhs);
  }

  void CVTermListInterface::replaceCVTerms(std::map<String, std::vector<CVTerm> > & cv_terms)
  {
    createIfNotExists_();
    cvt_ptr_->replaceCVTerms(cv_terms);
  }

  void CVTermListInterface::createIfNotExists_()
  {
    if (!cvt_ptr_) 
    {
      cvt_ptr_ = new CVTermList();
    }
  }

  bool CVTermListInterface::empty() const
  {
    return (cvt_ptr_ == nullptr || cvt_ptr_->empty());
  }

  void CVTermListInterface::setCVTerms(const std::vector<CVTerm>& terms)
  {
    createIfNotExists_();
    cvt_ptr_->setCVTerms(terms);
  }

  void CVTermListInterface::replaceCVTerm(const CVTerm& cv_term)
  {
    createIfNotExists_();
    cvt_ptr_->replaceCVTerm(cv_term);
  }

  void CVTermListInterface::replaceCVTerms(const std::vector<CVTerm>& cv_terms, const String& accession)
  {
    createIfNotExists_();
    cvt_ptr_->replaceCVTerms(cv_terms, accession);
  }

  void CVTermListInterface::replaceCVTerms(const std::map<String, std::vector<CVTerm> >& cv_term_map)
  {
    createIfNotExists_();
    cvt_ptr_->replaceCVTerms(cv_term_map);
  }

  void CVTermListInterface::consumeCVTerms(const std::map<String, std::vector<CVTerm> >& cv_term_map)
  {
    createIfNotExists_();
    cvt_ptr_->consumeCVTerms(cv_term_map);
  }

  const std::map<String, std::vector<CVTerm> >& CVTermListInterface::getCVTerms() const
  {
    if (!cvt_ptr_)
    {
      return empty_cvterm_map;
    }
    else
    {
      return cvt_ptr_->getCVTerms();
    }
  }

  /// adds a CV term
  void CVTermListInterface::addCVTerm(const CVTerm& term)
  {
    createIfNotExists_();
    cvt_ptr_->addCVTerm(term);
  }

  /// checks whether the term has a value
  bool CVTermListInterface::hasCVTerm(const String& accession) const
  {
    if (!cvt_ptr_)
    {
      return false;
    }
    else
    {
      return cvt_ptr_->hasCVTerm(accession);
    }
  }

}

