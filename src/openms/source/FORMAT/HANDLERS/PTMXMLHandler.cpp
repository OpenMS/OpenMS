// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

using namespace std;
using namespace xercesc;

namespace OpenMS::Internal
{

    PTMXMLHandler::PTMXMLHandler(map<std::string, pair<std::string, std::string> > & ptm_informations, const std::string & filename) :
      XMLHandler(filename, ""),
      ptm_informations_(ptm_informations)
    {
    }

    PTMXMLHandler::~PTMXMLHandler()
    = default;

    void PTMXMLHandler::writeTo(std::ostream & os)
    {
      os << "<PTMs>" << "\n";
      for (map<std::string, pair<std::string, std::string> >::const_iterator ptm_i = ptm_informations_.begin(); ptm_i != ptm_informations_.end(); ++ptm_i)
      {
        os << "\t<PTM>" << "\n";
        os << "\t\t<name>" << ptm_i->first << "</name>" << "\n";             // see header
        os << "\t\t<composition>" << ptm_i->second.first << "</composition>" << "\n";
        os << "\t\t<possible_amino_acids>" << ptm_i->second.second << "</possible_amino_acids>" << "\n";
        os << "\t</PTM>" << "\n";
      }
      os << "</PTMs>" << "\n";
    }

    void PTMXMLHandler::onStartElement(const char16_t * const qname, const XMLAttributes & /*attributes*/)
    {
      tag_ =StringUtils::trimmed(std::string(sm_.convert(qname)));
      open_tag_ = true;
    }

    void PTMXMLHandler::onEndElement(const char16_t* /*qname*/)
    {
//          tag_ =StringUtils::trimmed(StringUtils::toStr(sm_.convert(qname)));
      tag_ = "";
      open_tag_ = false;
    }

    void PTMXMLHandler::onCharacters(const char16_t* chars, Size /*length*/)
    {
      if (open_tag_)
      {
        if (tag_ == "name")
        {
          name_ =StringUtils::trimmed(std::string(sm_.convert(chars)));
        }
        else if (tag_ == "composition")
        {
          composition_ =StringUtils::trimmed(std::string(sm_.convert(chars)));
        }
        else if (tag_ == "possible_amino_acids")
        {
          ptm_informations_[name_] = make_pair(composition_,StringUtils::trimmed(std::string(sm_.convert(chars))));
        }
      }
    }
} // namespace OpenMS // namespace Internal
