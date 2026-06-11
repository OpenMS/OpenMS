// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>
#include <utility>

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

  DigestionEnzyme::DigestionEnzyme(const std::string& name,
                                   const std::string& cleavage_regex,
                                   const std::set<std::string>& synonyms,
                                   std::string regex_description) :
    name_(name),
    cleavage_regex_(cleavage_regex),
    synonyms_(synonyms),
    regex_description_(std::move(regex_description))
  {
  }

  DigestionEnzyme::DigestionEnzyme(const std::string& name,
                                   std::string cut_before,
                                   const std::string& nocut_after,
                                   std::string sense,
                                   const std::set<std::string>& synonyms,
                                   std::string regex_description) :
      name_(name),
      synonyms_(synonyms),
      regex_description_(std::move(regex_description))
  {
    //TODO check if all letters are A-Z?
    if (cut_before.empty())
    {
      //Maybe assertion?
      throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "No cleavage position given when trying to construct a DigestionEnzyme.");
    }
    else if (!StringUtils::hasSuffix(cut_before, "X"))
    {
      //TODO think about this
      cut_before = cut_before + "X";
    }
    cleavage_regex_ = "";
    if (StringUtils::toLower(sense) == "c")
    {
      cleavage_regex_ += "(?<=[" + cut_before + "]";
      if (!nocut_after.empty())
      {
        cleavage_regex_ += "(?!" + nocut_after + "])";
      }
    }
    else if (StringUtils::toLower(sense) == "n")
    {
      if (!nocut_after.empty())
      {
        cleavage_regex_ += "(?<![" + nocut_after + "])";
      }
      cleavage_regex_ += "(?=[" + cut_before + "]";
    }
    else
    {
      throw Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Cannot infer cleavage sense when constructing DigestionEnzyme. Has to be N or C.");
    }
  }

  DigestionEnzyme::~DigestionEnzyme() = default;

  void DigestionEnzyme::setName(const std::string& name)
  {
    name_ = name;
  }

  const std::string& DigestionEnzyme::getName() const
  {
    return name_;
  }

  void DigestionEnzyme::setSynonyms(const set<std::string>& synonyms)
  {
    synonyms_ = synonyms;
  }

  void DigestionEnzyme::addSynonym(const std::string& synonym)
  {
    synonyms_.insert(synonym);
  }

  const set<std::string>& DigestionEnzyme::getSynonyms() const
  {
    return synonyms_;
  }

  void DigestionEnzyme::setRegEx(const std::string& cleavage_regex)
  {
    cleavage_regex_ = cleavage_regex;
  }

  const std::string& DigestionEnzyme::getRegEx() const
  {
    return cleavage_regex_;
  }

  void DigestionEnzyme::setRegExDescription(const std::string& value)
  {
    regex_description_ = value;
  }

  const std::string& DigestionEnzyme::getRegExDescription() const
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

  bool DigestionEnzyme::operator==(const std::string& cleavage_regex) const
  {
    return cleavage_regex_ == cleavage_regex;
  }

  bool DigestionEnzyme::operator!=(const std::string& cleavage_regex) const
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

  bool DigestionEnzyme::setValueFromFile(const std::string& key, const std::string& value)
  {
    if (StringUtils::hasSuffix(key, ":Name"))
    {
      setName(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":RegEx"))
    {
      setRegEx(value);
      return true;
    }
    if (StringUtils::hasSuffix(key, ":RegExDescription"))
    {
      setRegExDescription(value);
      return true;
    }
    if (StringUtils::hasSubstring(key, ":Synonyms:"))
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

