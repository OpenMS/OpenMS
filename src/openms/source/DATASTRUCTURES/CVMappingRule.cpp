// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>

#include <OpenMS/DATASTRUCTURES/String.h>

using namespace std;

namespace OpenMS
{
  // CV mapping rule implementation
  CVMappingRule::CVMappingRule() :
    requirement_level_(CVMappingRule::MUST),
    combinations_logic_(CVMappingRule::OR)
  {
  }

  CVMappingRule::CVMappingRule(const CVMappingRule & rhs) :
    identifier_(rhs.identifier_),
    element_path_(rhs.element_path_),
    requirement_level_(rhs.requirement_level_),
    scope_path_(rhs.scope_path_),
    combinations_logic_(rhs.combinations_logic_),
    cv_terms_(rhs.cv_terms_)
  {
  }

  CVMappingRule::~CVMappingRule()
  {
  }

  CVMappingRule & CVMappingRule::operator=(const CVMappingRule & rhs)
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

  bool CVMappingRule::operator==(const CVMappingRule & rhs) const
  {
    return identifier_ == rhs.identifier_ &&
           element_path_ == rhs.element_path_ &&
           requirement_level_ == rhs.requirement_level_ &&
           scope_path_ == rhs.scope_path_ &&
           combinations_logic_ == rhs.combinations_logic_ &&
           cv_terms_ == rhs.cv_terms_;
  }

  bool CVMappingRule::operator!=(const CVMappingRule & rhs) const
  {
    return !(*this == rhs);
  }

  void CVMappingRule::setIdentifier(const String & identifier)
  {
    identifier_ = identifier;
  }

  const String & CVMappingRule::getIdentifier() const
  {
    return identifier_;
  }

  void CVMappingRule::setElementPath(const String & element_path)
  {
    element_path_ = element_path;
  }

  const String & CVMappingRule::getElementPath() const
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

  void CVMappingRule::setScopePath(const String & path)
  {
    scope_path_ = path;
  }

  const String & CVMappingRule::getScopePath() const
  {
    return scope_path_;
  }

  void CVMappingRule::setCVTerms(const vector<CVMappingTerm> & cv_terms)
  {
    cv_terms_ = cv_terms;
  }

  const vector<CVMappingTerm> & CVMappingRule::getCVTerms() const
  {
    return cv_terms_;
  }

  void CVMappingRule::addCVTerm(const CVMappingTerm & cv_term)
  {
    cv_terms_.push_back(cv_term);
  }

} // namespace OpenMS
