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

#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iostream>


namespace OpenMS {
class CVMappingRule;
}  // namespace OpenMS

using namespace std;

namespace OpenMS
{
  CVMappings::CVMappings()
  {
  }

  CVMappings::CVMappings(const CVMappings& rhs) :
    mapping_rules_(rhs.mapping_rules_),
    cv_references_(rhs.cv_references_),
    cv_references_vector_(rhs.cv_references_vector_)
  {
  }

  CVMappings::~CVMappings()
  {
  }

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
