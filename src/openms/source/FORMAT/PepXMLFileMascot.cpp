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
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/PepXMLFileMascot.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS
{

  PepXMLFileMascot::PepXMLFileMascot() :
    XMLHandler("", "1.8"),
    XMLFile("/SCHEMAS/PepXML_1_8.xsd", "1.8"),
    peptides_(nullptr)
  {

  }

  void PepXMLFileMascot::load(const String & filename, map<String, vector<AASequence> > & peptides)
  {
    //Filename for error messages in XMLHandler
    file_ = filename;

    peptides.clear();

    peptides_ = &peptides;

    parse_(filename, this);

    //reset members
    actual_title_ = "";
    actual_sequence_ = "";
    actual_modifications_ = vector<pair<String, UInt> >();
    peptides_ = nullptr;
    variable_modifications_ = vector<pair<String, double> >();
    fixed_modifications_ = vector<String>();
  }

  void PepXMLFileMascot::matchModification_(double mass, String & modification_description)
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
    String element = sm_.convert(qname);

    //cout << "Start: " << element << "\n";

    //SEARCH PARAMETERS
    if (element == "aminoacid_modification")
    {
      String temp_string = attributeAsString_(attributes, "variable");
      if (temp_string == "Y")
      {
        variable_modifications_.push_back(make_pair(attributeAsString_(attributes, "description"),
                                                    attributeAsDouble_(attributes, "mass")));
      }
      else
      {
        fixed_modifications_.push_back(attributeAsString_(attributes, "description"));
      }
    }

    // <terminal_modification terminus="n" massdiff="+108.05" mass="109.06" variable="N" protein_terminus="" description="dNIC (N-term)"/>
    if (element == "terminal_modification")
    {
      String temp_string = attributeAsString_(attributes, "variable");
      if (temp_string == "Y")
      {
        variable_modifications_.push_back(make_pair(attributeAsString_(attributes, "description"),
                                                    attributeAsDouble_(attributes, "mass")));

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
      String temp_description = "";
      UInt modification_position = attributeAsInt_(attributes, "position");
      double modification_mass = attributeAsDouble_(attributes, "mass");

      matchModification_(modification_mass, temp_description);

      // the modification position is 1-based
      actual_modifications_.push_back(make_pair(temp_description, modification_position));
    }
  }

  void PepXMLFileMascot::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
  {
    String element = sm_.convert(qname);

    ///SEARCH PARAMETERS
    if (element == "search_hit")
    {
      AASequence temp_aa_sequence = AASequence::fromString(actual_sequence_);

      // modification position is 1-based
      for (vector<pair<String, UInt> >::const_iterator it = actual_modifications_.begin(); it != actual_modifications_.end(); ++it)
      {
        // e.g. Carboxymethyl (C)
        vector<String> mod_split;
        it->first.split(' ', mod_split);
        if (it->first.hasSubstring("C-term"))
        {
          temp_aa_sequence.setCTerminalModification(it->first);
        }
        else if (it->first.hasSubstring("N-term"))
        {
          temp_aa_sequence.setNTerminalModification(it->first);
        }
        if (mod_split.size() == 2)
        {
          temp_aa_sequence.setModification(it->second - 1, mod_split[0]);
        }
        else
        {
          error(LOAD, String("Cannot parse modification '") + it->first + "@" + it->second + "'");
        }
      }

      // fixed modifications
      for (vector<String>::const_iterator it = fixed_modifications_.begin(); it != fixed_modifications_.end(); ++it)
      {
        // e.g. Carboxymethyl (C)
        vector<String> mod_split;
        it->split(' ', mod_split);
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
              String origin = mod_split[1];
              origin.remove(')');
              origin.remove('(');
              for (Size i = 0; i != temp_aa_sequence.size(); ++i)
              {
                // best way we can check; because origin can be e.g. (STY)
                if (origin.hasSubstring(temp_aa_sequence[i].getOneLetterCode()))
                {
                  temp_aa_sequence.setModification(i, mod_split[0]);
                }
              }
            }
          }
        }
        else
        {
          error(LOAD, String("Cannot parse fixed modification '") + *it + "'");
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
