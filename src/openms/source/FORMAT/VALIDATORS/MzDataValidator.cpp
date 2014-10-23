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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/VALIDATORS/MzDataValidator.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    MzDataValidator::MzDataValidator(const CVMappings & mapping, const ControlledVocabulary & cv) :
      SemanticValidator(mapping, cv)
    {
      setCheckUnits(true);
    }

    MzDataValidator::~MzDataValidator()
    {
    }

    void MzDataValidator::handleTerm_(const String & path, const CVTerm & parsed_term)
    {
      //check if the term is allowed in this element
      //and if there is a mapping rule for this element
      //Also store fulfilled rule term counts - this count is used to check of the MUST/MAY and AND/OR/XOR is fulfilled
      bool allowed = false;
      bool rule_found = false;
      vector<CVMappingRule> & rules = rules_[path];
      for (Size r = 0; r < rules.size(); ++r)  //go thru all rules
      {
        rule_found = true;
        for (Size t = 0; t < rules[r].getCVTerms().size(); ++t)     //go thru all terms
        {
          const CVMappingTerm & term = rules[r].getCVTerms()[t];
          if (term.getUseTerm() && term.getAccession() == parsed_term.accession)         //check if the term itself is allowed
          {
            allowed = true;
            fulfilled_[path][rules[r].getIdentifier()][term.getAccession()]++;
            break;
          }
          if (term.getAllowChildren())           //check if the term's children are allowed
          {
            set<String> child_terms;
            cv_.getAllChildTerms(child_terms, term.getAccession());
            for (set<String>::const_iterator it = child_terms.begin(); it != child_terms.end(); ++it)
            {
              if (*it == parsed_term.accession)
              {
                allowed = true;
                fulfilled_[path][rules[r].getIdentifier()][term.getAccession()]++;
                break;
              }
            }
          }
        }
      }

      // check units
      if (check_units_ && cv_.exists(parsed_term.accession))
      {
        ControlledVocabulary::CVTerm term = cv_.getTerm(parsed_term.accession);
        // check if the cv term has units
        if (!term.units.empty())
        {
          if (!parsed_term.has_unit_accession)
          {
            errors_.push_back(String("CV term must have a unit: " + parsed_term.accession + " - " + parsed_term.name));
          }
          else
          {
            // check if the accession is ok
            if (cv_.exists(parsed_term.unit_accession))
            {
              // check whether this unit is allowed within the cv term
              if (term.units.find(parsed_term.unit_accession) == term.units.end())
              {
                // last chance, a child term of the units was used
                set<String> child_terms;

                bool found_unit(false);
                for (set<String>::const_iterator it = term.units.begin(); it != term.units.end(); ++it)
                {
                  cv_.getAllChildTerms(child_terms, *it);
                  if (child_terms.find(parsed_term.unit_accession) != child_terms.end())
                  {
                    found_unit = true;
                    break;
                  }
                }

                if (!found_unit)
                {
                  errors_.push_back(String("Unit CV term not allowed: " + parsed_term.unit_accession + " - " + parsed_term.unit_name + " of term " + parsed_term.accession + " - " + parsed_term.name));
                }
              }
            }
            else
            {
              errors_.push_back(String("Unit CV term not found: " + parsed_term.unit_accession + " - " + parsed_term.unit_name + " of term " + parsed_term.accession + " - " + parsed_term.name));
            }
          }
        }
        else
        {
          // check whether unit was used
          if (parsed_term.has_unit_accession || parsed_term.has_unit_name)
          {
            warnings_.push_back(String("Unit CV term used, but not allowed: " + parsed_term.unit_accession + " - " + parsed_term.unit_name + " of term " + parsed_term.accession + " - " + parsed_term.name));
          }
        }
      }

      if (!rule_found)       //No rule found
      {
        warnings_.push_back(String("No mapping rule found for element '") + getPath_(1) + "'");
      }
      else if (!allowed)      //if rule found and not allowed
      {
        errors_.push_back(String("CV term used in invalid element: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
      }

      //check if term accession and term name match
      if (cv_.exists(parsed_term.accession))
      {
        String parsed_name = parsed_term.name;
        parsed_name.trim();
        String correct_name = cv_.getTerm(parsed_term.accession).name;
        correct_name.trim();

        //be a bit more soft: ignore upper-lower case
        parsed_name.toLower();
        correct_name.toLower();

        //be a bit more soft: ignore spaces
        parsed_name.removeWhitespaces();
        correct_name.removeWhitespaces();

        if (parsed_name != correct_name)
        {
          errors_.push_back(String("Name of CV term not correct: '") + parsed_term.accession + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
      }
    }

  }   // namespace Internal
} // namespace OpenMS
