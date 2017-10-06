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
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <iostream>

using namespace std;

namespace OpenMS
{
  DigestionEnzyme::DigestionEnzyme() :
    name_("unknown_enzyme"),
    cleavage_regex_(""),
    synonyms_(),
    regex_description_("")
  {
  }

  DigestionEnzyme::DigestionEnzyme(const DigestionEnzyme& enzyme) :
    name_(enzyme.name_),
    cleavage_regex_(enzyme.cleavage_regex_),
    synonyms_(enzyme.synonyms_),
    regex_description_(enzyme.regex_description_)
  {
  }

  DigestionEnzyme::DigestionEnzyme(const String& name,
                                   const String& cleavage_regex,
                                   const std::set<String>& synonyms,
                                   String regex_description) :
    name_(name),
    cleavage_regex_(cleavage_regex),
    synonyms_(synonyms),
    regex_description_(regex_description)
  {
  }

  DigestionEnzyme::~DigestionEnzyme()
  {
  }

  DigestionEnzyme& DigestionEnzyme::operator=(const DigestionEnzyme& enzyme)
  {
    if (this != &enzyme)
    {
      name_ = enzyme.name_;
      cleavage_regex_ = enzyme.cleavage_regex_;
      synonyms_ = enzyme.synonyms_;
      regex_description_ = enzyme.regex_description_;
    }
    return *this;
  }

  void DigestionEnzyme::setName(const String& name)
  {
    name_ = name;
  }

  String DigestionEnzyme::getName() const
  {
    return name_;
  }

  void DigestionEnzyme::setSynonyms(const set<String>& synonyms)
  {
    synonyms_ = synonyms;
  }

  void DigestionEnzyme::addSynonym(const String& synonym)
  {
    synonyms_.insert(synonym);
  }

  const set<String>& DigestionEnzyme::getSynonyms() const
  {
    return synonyms_;
  }

  void DigestionEnzyme::setRegEx(const String& cleavage_regex)
  {
    cleavage_regex_ = cleavage_regex;
  }

  String DigestionEnzyme::getRegEx() const
  {
    return cleavage_regex_;
  }

  void DigestionEnzyme::setRegExDescription(const String& value)
  {
    regex_description_ = value;
  }

  String DigestionEnzyme::getRegExDescription() const
  {
    return regex_description_;
  }

  bool DigestionEnzyme::operator==(const DigestionEnzyme& enzyme) const
  {
    return name_ == enzyme.name_ &&
           synonyms_ == enzyme.synonyms_ &&
           cleavage_regex_ == enzyme.cleavage_regex_ &&
           regex_description_ == enzyme.regex_description_;
  }

  bool DigestionEnzyme::operator==(const String& cleavage_regex) const
  {
    return cleavage_regex_ == cleavage_regex;
  }

  bool DigestionEnzyme::operator!=(const String& cleavage_regex) const
  {
    return cleavage_regex_ != cleavage_regex;
  }

  bool DigestionEnzyme::operator!=(const DigestionEnzyme& enzyme) const
  {
    return !(*this == enzyme);
  }

  bool DigestionEnzyme::operator<(const DigestionEnzyme& enzyme) const
  {
    return this->getName() < enzyme.getName();
  }

  bool DigestionEnzyme::setValueFromFile(const String& key, const String& value)
  {
    if (key.hasSuffix(":Name"))
    {
      setName(value);
      return true;
    }
    if (key.hasSuffix(":RegEx"))
    {
      setRegEx(value);
      return true;
    }
    if (key.hasSuffix(":RegExDescription"))
    {
      setRegExDescription(value);
      return true;
    }
    if (key.hasSubstring(":Synonyms:"))
    {
      addSynonym(value);
      return true;
    }
    return false;
  }

  ostream& operator<<(ostream& os, const DigestionEnzyme& enzyme)
  {
    os << "digestion enzyme:" << enzyme.name_ << " (cleavage: "
       << enzyme.cleavage_regex_ << " - " << enzyme.regex_description_ << ")";
    return os;
  }

}
