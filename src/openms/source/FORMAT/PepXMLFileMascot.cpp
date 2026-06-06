// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PepXMLFileMascot.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS
{

  PepXMLFileMascot::PepXMLFileMascot() :
    XMLHandler("", "1.8"),
    XMLFile("/SCHEMAS/PepXML_1_8.xsd", "1.8"),
    peptides_(nullptr)
  {

  }

  void PepXMLFileMascot::load(const std::string & filename, map<std::string, vector<AASequence> > & peptides)
  {
    //Filename for error messages in XMLHandler
    file_ = filename;

    peptides.clear();

    peptides_ = &peptides;

    parse_(filename, this);

    //reset members
    actual_title_ = "";
    actual_sequence_ = "";
    actual_modifications_ = vector<pair<std::string, UInt> >();
    peptides_ = nullptr;
    variable_modifications_ = vector<pair<std::string, double> >();
    fixed_modifications_ = vector<std::string>();
  }

  void PepXMLFileMascot::matchModification_(double mass, std::string & modification_description)
  {
    UInt i = 0;
    bool found = false;

    while (i < variable_modifications_.size() && !found)
    {
      double difference = variable_modifications_[i].second - mass;
      if (difference < 0)
      {
        difference *= -1;
      }
      if (difference < 0.001)
      {
        modification_description = variable_modifications_[i].first;
        found = true;
      }
      ++i;
    }
  }

  void PepXMLFileMascot::startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes)
  {
    std::string element = sm_.convert(qname);

    //cout << "Start: " << element << "\n";

    //SEARCH PARAMETERS
    if (element == "aminoacid_modification")
    {
      std::string temp_string = attributeAsString_(attributes, "variable");
      if (temp_string == "Y")
      {
        variable_modifications_.emplace_back(attributeAsString_(attributes, "description"),
                                                    attributeAsDouble_(attributes, "mass"));
      }
      else
      {
        fixed_modifications_.push_back(attributeAsString_(attributes, "description"));
      }
    }

    // <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
    if (element == "terminal_modification")
    {
      std::string temp_string = attributeAsString_(attributes, "variable");
      if (temp_string == "Y")
      {
        variable_modifications_.emplace_back(attributeAsString_(attributes, "description"),
                                                    attributeAsDouble_(attributes, "mass"));

      }
      else
      {
        fixed_modifications_.push_back(attributeAsString_(attributes, "description"));
      }
    }
    //PEPTIDES
    else if (element == "spectrum_query")
    {
      actual_title_ = attributeAsString_(attributes, "spectrum");
    }
    else if (element == "search_hit")
    {
      actual_sequence_ = attributeAsString_(attributes, "peptide");
    }
    else if (element == "mod_aminoacid_mass")
    {
      std::string temp_description = "";
      UInt modification_position = attributeAsInt_(attributes, "position");
      double modification_mass = attributeAsDouble_(attributes, "mass");

      matchModification_(modification_mass, temp_description);

      // the modification position is 1-based
      actual_modifications_.emplace_back(temp_description, modification_position);
    }
  }

  void PepXMLFileMascot::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
  {
    std::string element = sm_.convert(qname);

    ///SEARCH PARAMETERS
    if (element == "search_hit")
    {
      AASequence temp_aa_sequence = AASequence::fromString(actual_sequence_);

      // modification position is 1-based
      for (vector<pair<std::string, UInt> >::const_iterator it = actual_modifications_.begin(); it != actual_modifications_.end(); ++it)
      {
        // e.g. Carboxymethyl (C)
        vector<std::string> mod_split;
        StringUtils::split(it->first, ' ', mod_split);
        if (StringUtils::hasSubstring(it->first, "C-term"))
        {
          temp_aa_sequence.setCTerminalModification(it->first);
        }
        else if (StringUtils::hasSubstring(it->first, "N-term"))
        {
          temp_aa_sequence.setNTerminalModification(it->first);
        }
        if (mod_split.size() == 2)
        {
          temp_aa_sequence.setModification(it->second - 1, mod_split[0]);
        }
        else
        {
          error(LOAD,StringUtils::toStr("Cannot parse modification '") + it->first + "@" + it->second + "'");
        }
      }

      // fixed modifications
      for (const std::string& it : fixed_modifications_)
      {
        // e.g. Carboxymethyl (C)
        vector<std::string> mod_split;
        StringUtils::split(it, ' ', mod_split);
        if (mod_split.size() == 2)
        {
          if (mod_split[1] == "(C-term)")
          {
            temp_aa_sequence.setCTerminalModification(mod_split[0]);
          }
          else
          {
            if (mod_split[1] == "(N-term)")
            {
              temp_aa_sequence.setNTerminalModification(mod_split[0]);
            }
            else
            {
              std::string origin = mod_split[1];
              StringUtils::remove(origin, ')');
              StringUtils::remove(origin, '(');
              for (Size i = 0; i != temp_aa_sequence.size(); ++i)
              {
                // best way we can check; because origin can be e.g. (STY)
                if (StringUtils::hasSubstring(origin, temp_aa_sequence[i].getOneLetterCode()))
                {
                  temp_aa_sequence.setModification(i, mod_split[0]);
                }
              }
            }
          }
        }
        else
        {
          error(LOAD,StringUtils::toStr("Cannot parse fixed modification '") + it + "'");
        }
      }

      actual_aa_sequences_.push_back(temp_aa_sequence);

      actual_modifications_.clear();
    }
    else if (element == "spectrum_query")
    {
      peptides_->insert(make_pair(actual_title_, actual_aa_sequences_));
      actual_aa_sequences_.clear();
    }
  }

} // namespace OpenMS
