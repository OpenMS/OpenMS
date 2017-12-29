// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTermListInterface.h>

#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/Helpers.h>

#include <iostream>

namespace OpenMS
{

  static const Map<String, std::vector<CVTerm> > empty_cvterm_map = Map<String, std::vector<CVTerm> >();

  CVTermListInterface::CVTermListInterface() :
      MetaInfoInterface(),
      cvt_ptr_(nullptr)
    {}


  CVTermListInterface::CVTermListInterface(const CVTermListInterface & rhs) :
    MetaInfoInterface(rhs),
    cvt_ptr_(nullptr)
  {
    MetaInfoInterface::operator=(rhs);

    if (rhs.cvt_ptr_ != nullptr)
    {
      cvt_ptr_ = new CVTermList(*rhs.cvt_ptr_);
    }
  }

    // Destructor (non virtual)
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

    bool CVTermListInterface::operator==(const CVTermListInterface& rhs) const
    {
      return MetaInfoInterface::operator==(rhs) &&
             Helpers::cmpPtrSafe<CVTermList*>(cvt_ptr_, rhs.cvt_ptr_);
    }

    bool CVTermListInterface::operator!=(const CVTermListInterface& rhs) const
    {
      return !(*this == rhs);
    }

    void CVTermListInterface::replaceCVTerms(Map<String, std::vector<CVTerm> > & cv_terms)
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

  void CVTermListInterface::replaceCVTerms(const Map<String, std::vector<CVTerm> >& cv_term_map)
  {
    createIfNotExists_();
    cvt_ptr_->replaceCVTerms(cv_term_map);
  }

  void CVTermListInterface::consumeCVTerms(const Map<String, std::vector<CVTerm> >& cv_term_map)
  {
    createIfNotExists_();
    cvt_ptr_->consumeCVTerms(cv_term_map);
  }

  const Map<String, std::vector<CVTerm> >& CVTermListInterface::getCVTerms() const
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
