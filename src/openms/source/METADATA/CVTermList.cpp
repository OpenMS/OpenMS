// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/CVTermList.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace std;

namespace OpenMS
{

  CVTermList::~CVTermList()
  {
  }

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
    for (vector<CVTerm>::const_iterator it = cv_terms.begin(); it != cv_terms.end(); ++it)
    {
      addCVTerm(*it);
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

  void CVTermList::replaceCVTerms(const Map<String, vector<CVTerm> >& cv_term_map)
  {
    cv_terms_ = cv_term_map;
  }

  void CVTermList::consumeCVTerms(const Map<String, vector<CVTerm> >& cv_term_map)
  {
    for (std::map<String, std::vector<CVTerm> >::const_iterator it = cv_term_map.begin(); it != cv_term_map.end(); ++it)
    {
      cv_terms_[it->first].insert(cv_terms_[it->first].end(), it->second.begin(), it->second.end());
    }
  }

  const Map<String, vector<CVTerm> >& CVTermList::getCVTerms() const
  {
    return cv_terms_;
  }

  bool CVTermList::hasCVTerm(const String& accession) const
  {
    return cv_terms_.has(accession);
  }

  bool CVTermList::operator==(const CVTermList& cv_term_list) const
  {
    return MetaInfoInterface::operator==(cv_term_list) && cv_terms_.equals(cv_term_list.cv_terms_);
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

