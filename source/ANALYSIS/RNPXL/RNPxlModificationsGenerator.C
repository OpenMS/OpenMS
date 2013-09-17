// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

using namespace std;

#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
namespace OpenMS
{

//static
bool RNPxlModificationsGenerator::notInSeq(String res_seq, String query)
{
  // special case: empty query is in every seq -> false
  if (query == "")
  {
    return false;
  }

  // test all k-mers with k=size of query
  for (Int l = 0; l <= (Int)res_seq.size() - (Int)query.size(); ++l)
  {
    String a = res_seq.substr(l, query.size());
    String b = query;

    sort(a.begin(), a.end());
    sort(b.begin(), b.end());

    if (a == b)
    {
      return false;
    }
  }
  return true;
}

//static
RNPxlModificationMassesResult RNPxlModificationsGenerator::initModificationMassesRNA(StringList target_nucleotides, StringList mappings, StringList restrictions, StringList modifications, String sequence_restriction, bool cysteine_adduct, Int max_length)
{
  String original_sequence_restriction = sequence_restriction;

  // 152 modification
  const String cysteine_adduct_string("C1H2N3O6");
  const EmpiricalFormula cysteine_adduct_formula(cysteine_adduct_string); // 152 modification

  RNPxlModificationMassesResult result;
  // read nucleotides and empirical formula of monophosphat version
  map<String, EmpiricalFormula> map_target_to_formula;
  for (Size i = 0; i != target_nucleotides.size(); ++i)
  {
    String s = target_nucleotides[i];
    vector<String> fields;
    s.split("=", fields);
    map_target_to_formula[fields[0]] = EmpiricalFormula(fields[1]);
  }

  // read mapping of source to target
  map<char, vector<char> > map_source_to_targets;
  for (Size i = 0; i != mappings.size(); ++i)
  {
    String s = mappings[i];
    vector<String> fields;
    s.split("->", fields);
    map_source_to_targets[fields[0][0]].push_back(fields[1][0]);
  }

  // extract source nucleotides based on mapping
  vector<char> source_nucleotides; // nucleotides as expected in the restriction sequence
  for (StringList::const_iterator sit = mappings.begin(); sit != mappings.end(); ++sit)
  {
    source_nucleotides.push_back((*sit)[0]);
  }

  if (sequence_restriction.empty())
  {
    vector<String> all_combinations;
    vector<String> actual_combinations;

    // add single source nucleotides to all_combinations
    for (Size i = 0; i != source_nucleotides.size(); ++i)
    {
      all_combinations.push_back(String(source_nucleotides[i]));
      actual_combinations.push_back(String(source_nucleotides[i]));
    }

    for (Int i = 1; i <= max_length - 1; ++i)
    {
      vector<String> new_combinations;
      for (Size n = 0; n != source_nucleotides.size(); ++n)
      {
        // grow actual_combinations/ all_combinations by one nucleotide
        for (Size c = 0; c != actual_combinations.size(); ++c)
        {
          new_combinations.push_back(source_nucleotides[n] + actual_combinations[c]);
          all_combinations.push_back(source_nucleotides[n] + actual_combinations[c]);
        }
      }
      actual_combinations = new_combinations;
    }

    for (Size i = 0; i != all_combinations.size(); ++i)
    {
      sequence_restriction += all_combinations[i];
    }
  }

  // read restrictions
  std::cout << "Min. count restrictions:" << endl;
  map<char, Size> map_target_to_mincount;
  for (Size i = 0; i != restrictions.size(); ++i)
  {
    String s = restrictions[i];
    vector<String> fields;
    s.split("=", fields);
    Size min_count = (Size)fields[1].toInt();
    if (min_count > 0)
    {
      map_target_to_mincount[fields[0][0]] = min_count;
      std::cout << "\t" << "min. count: " << fields[0][0] << "\t" << min_count << endl;
    }
  }

  //cout << "source sequence: " << sequence_restriction << endl;

  // erase trivial cases from mapping so only the combinatorial cases: 1 source -> n targets remain
  for (map<char, vector<char> >::iterator sit = map_source_to_targets.begin(); sit != map_source_to_targets.end(); )
  {
    char source = sit->first;
    char first_target = sit->second[0];

    if (sit->second.size() == 1 && source == first_target) // trivial case e.g. A->A... no substitution needed
    {
      map_source_to_targets.erase(sit++);
    }
    else if (sit->second.size() == 1 && source != first_target) // simple rename e.g. A->X... simply subsitute all in restriction sequence
    {
      sequence_restriction.substitute(source, first_target);
      map_source_to_targets.erase(sit++);
    }
    else // multiple targets
    {
      ++sit;
    }
  }

  // std::cout << "source sequence: " << sequence_restriction << endl;

  if (!map_source_to_targets.empty() && sequence_restriction.empty())
  {
    std::cout << "WARNING: no restriction on sequence but multiple target nucleotides specified. Will generate huge amount of sequences" << endl;
  }

  vector<vector<bool> > modifications_is_subtractive(modifications.size(), vector<bool>());
  vector<vector<EmpiricalFormula> > modification_formulas(modifications.size(), vector<EmpiricalFormula>());
  for (Size i = 0; i != modifications.size(); ++i)
  {
    // decompose string into additive and subtractive EmpiricalFormulas
    modifications[i].substitute("-", "#-");
    modifications[i].substitute("+", "#+");
    vector<String> ems;
    modifications[i].split("#", ems);

    for (Size j = 0; j != ems.size(); ++j)
    {
      if (ems[j].empty())
      {
        continue;
      }

      if (ems[j][0] == '-')
      {
        modifications_is_subtractive[i].push_back(true);
        ems[j].remove('-');
      }
      else if (ems[j][0] == '+')
      {
        modifications_is_subtractive[i].push_back(false);
        ems[j].remove('+');
      }
      else // no + still means additive
      {
        modifications_is_subtractive[i].push_back(false);
      }
      EmpiricalFormula ef(ems[j]);
      ef.setCharge(0);
      modification_formulas[i].push_back(ef);
    }

    std::cout << "Modification: " << endl;
    for (Size f = 0; f != modification_formulas[i].size(); ++f)
    {
      std::cout << "\t" << modification_formulas[i][f] << " subtractive: " << modifications_is_subtractive[i][f] << endl;
    }
  }

  // generate all target sequences by substituting each source nucleotide by their target nucleotide(s)
  StringList target_sequences;
  generateTargetSequences(sequence_restriction, 0, map_source_to_targets, target_sequences);
  std::cout << "target sequence(s):" << target_sequences.size() << endl;

  if (!original_sequence_restriction.empty())
  {
    for (Size i = 0; i != target_sequences.size(); ++i)
    {
      if (target_sequences[i].size() < 60)
      {
        std::cout << target_sequences[i] << endl;
      }
      else
      {
        std::cout << target_sequences[i].prefix(60) << "..."  << endl;
      }
    }
  }

  {
    vector<EmpiricalFormula> actual_combinations;
    for (map<String, EmpiricalFormula>::const_iterator mit = map_target_to_formula.begin(); mit != map_target_to_formula.end(); ++mit)
    {
      String target_nucleotide = mit->first;
      std::cout << "target nucleotide: " << target_nucleotide << endl;
      EmpiricalFormula target_nucleotide_formula = mit->second;
      for (Size m = 0; m != modification_formulas.size(); ++m)
      {
        EmpiricalFormula e(target_nucleotide_formula);
        String s(target_nucleotide);
        for (Size f = 0; f != modification_formulas[m].size(); ++f)
        {
          if (modifications_is_subtractive[m][f])
          {
            e = e - modification_formulas[m][f];
            s += "-" + modification_formulas[m][f].getString();
          }
          else
          {
            e = e + modification_formulas[m][f];
            s += "+" + modification_formulas[m][f].getString();
          }
        }
        actual_combinations.push_back(e);
        result.mod_combinations[actual_combinations.back().getString()].insert(s);

        std::cout << "\t" << "modifications: " << s << "\t\t" << e.getString() << endl;
      }
    }

    vector<EmpiricalFormula> all_combinations = actual_combinations;

    for (Int i = 0; i < max_length - 1; ++i)
    {
      vector<EmpiricalFormula> new_combinations;
      for (map<String, EmpiricalFormula>::const_iterator mit = map_target_to_formula.begin(); mit != map_target_to_formula.end(); ++mit)
      {
        String target_nucleotide = mit->first;
        EmpiricalFormula target_nucleotide_formula = mit->second;
        for (Size c = 0; c != actual_combinations.size(); ++c)
        {
          new_combinations.push_back(target_nucleotide_formula + actual_combinations[c] - EmpiricalFormula("H2O")); // -H2O because of condensation reaction
          all_combinations.push_back(target_nucleotide_formula + actual_combinations[c] - EmpiricalFormula("H2O")); // " "
          const set<String>& ambiguities = result.mod_combinations[actual_combinations[c].getString()];
          for (set<String>::const_iterator sit = ambiguities.begin(); sit != ambiguities.end(); ++sit)
          {
            result.mod_combinations[all_combinations.back().getString()].insert(target_nucleotide + *sit);
          }
          //cout << target_nucleotide + mod_combinations[actual_combinations[c].getString()]  << endl;
        }
      }
      actual_combinations = new_combinations;
    }

    //    std::cout << all_combinations.size() << endl;
    for (Size i = 0; i != all_combinations.size(); ++i)
    {
      //      std::cout << all_combinations[i].getString() << endl;
      result.mod_masses[all_combinations[i].getString()] = all_combinations[i].getMonoWeight();
    }
  }

  std::cout << "Filtering on restrictions... " << endl;
  // filtering on restrictions
  std::vector<pair<String, String> > violates_restriction; // elemental composition, nucleotide style formula
  for (Map<String, DoubleReal>::ConstIterator mit = result.mod_masses.begin(); mit != result.mod_masses.end(); ++mit)
  {
    // remove additive or subtractive modifications from string as these are not used in string comparison
    const set<String>& ambiguities = result.mod_combinations[mit->first];
    for (set<String>::const_iterator sit = ambiguities.begin(); sit != ambiguities.end(); ++sit)
    {
      String nucleotide_style_formula = *sit;
      Size p1 = nucleotide_style_formula.find('-');
      Size p2 = nucleotide_style_formula.find('+');
      Size p = min(p1, p2);
      if (p != String::npos)
      {
        nucleotide_style_formula = nucleotide_style_formula.prefix(p);
      }
      //  std::cout << "(" << nucleotide_style_formula << ")" << endl;
      // perform string comparison: if nucleotide seuquence doesnt occur in any permutation in res_seq mark the corresponding empirical formula for deletion
      bool restriction_violated = false;

      // for each min. count restriction on a target nucleotide...
      for (map<char, Size>::const_iterator minit = map_target_to_mincount.begin(); minit != map_target_to_mincount.end(); ++minit)
      {
        //    std::cout << nucleotide_style_formula <<  " current target: " << minit->first << " ";
        Size occurances = (Size) std::count(nucleotide_style_formula.begin(), nucleotide_style_formula.end(), minit->first);
        //    std::cout << occurances << endl;
        if (occurances < minit->second)
        {
          restriction_violated = true;
        }
      }

      // check if contained in at least one of the target sequences
      bool containment_violated = false;
      Size violation_count = 0;
      for (StringList::const_iterator tsit = target_sequences.begin(); tsit != target_sequences.end(); ++tsit)
      {
        String current_target_seq = *tsit;
        if (notInSeq(current_target_seq, nucleotide_style_formula))
        {
          ++violation_count;
        }
      }

      if (violation_count == target_sequences.size())
      {
        containment_violated = true;
      }

      if (containment_violated || restriction_violated)
      {
        violates_restriction.push_back(make_pair(mit->first, *sit)); // chemical formula, nucleotide style formula pair violates restrictions
      }
    }
  }

  for (size_t i = 0; i != violates_restriction.size(); ++i)
  {
    const String& chemical_formula = violates_restriction[i].first;
    result.mod_combinations[chemical_formula].erase(violates_restriction[i].second);
    //cout << "filtered sequence: " << chemical_formula << "\t" << violates_restriction[i].first << endl;
  }

  // standard associative-container erase idiom
  for (map<String, set<String> >::iterator mcit = result.mod_combinations.begin(); mcit != result.mod_combinations.end(); )
  {
    if (mcit->second.empty())
    {
      //cout << "filtered sequence: " << mcit->first << endl;
      result.mod_masses.erase(mcit->first); // remove from mod masses
      result.mod_combinations.erase(mcit++); // don't change precedence !
    }
    else
    {
      ++mcit;   // don't change precedence !
    }
  }

  if (cysteine_adduct)
  {
    result.mod_masses[cysteine_adduct_formula.getString()] = cysteine_adduct_formula.getMonoWeight();
    result.mod_combinations[cysteine_adduct_formula.getString()].insert(cysteine_adduct_string);
  }

  // output index  -> empirical formula -> (ambigous) nucleotide formulas
  // nucleotide formulas which only differ in nucleotide ordering are only printed once
  DoubleReal pseudo_rt = 1;
  for (Map<String, DoubleReal>::ConstIterator mit = result.mod_masses.begin(); mit != result.mod_masses.end(); ++mit)
  {
    result.mod_formula_idx[pseudo_rt] = mit->first;

    if (cysteine_adduct && mit->first == cysteine_adduct_formula.getString())
    {
      std::cout << pseudo_rt++ << " " << mit->first << " " << mit->second << " ( cysteine adduct )" << endl;
      continue;
    }

    std::cout << pseudo_rt++ << " " << mit->first << " " << mit->second << " ( ";

    const set<String>& ambiguities = result.mod_combinations[mit->first];
    set<String> printed;
    for (set<String>::const_iterator sit = ambiguities.begin(); sit != ambiguities.end(); ++sit)
    {
      String nucleotide_style_formula = *sit;
      Size p1 = nucleotide_style_formula.find('-');
      Size p2 = nucleotide_style_formula.find('+');
      Size p = min(p1, p2);
      // sort nucleotides up to beginning of modification (first '+' or '-')
      if (p != String::npos)
      {
        std::sort(nucleotide_style_formula.begin(), nucleotide_style_formula.begin() + p);
      }
      else
      {
        std::sort(nucleotide_style_formula.begin(), nucleotide_style_formula.end());
      }

      // only print ambigous sequences once
      if (printed.find(nucleotide_style_formula) == printed.end())
      {
        std::cout << nucleotide_style_formula;
        std::cout << " ";
        printed.insert(nucleotide_style_formula);
      }
    }

    std::cout << ")" << endl;
  }
  std::cout << "Finished generation of modification masses." << endl;
  return result;
}

//static
void  RNPxlModificationsGenerator::generateTargetSequences(const String& res_seq, Size pos, const map<char, vector<char> >& map_source2target, StringList& target_sequences)
{
  typedef map<char, vector<char> >::const_iterator TConstMapIterator;

  while (pos < res_seq.size())
  {
    // check if current character is in source 2 target map
    TConstMapIterator target_iterator = map_source2target.find(res_seq[pos]);
    if (target_iterator == map_source2target.end())
    {
      ++pos;
    }
    else // yes?
    {
      const vector<char>& targets = target_iterator->second;
      for (Size i = 0; i != targets.size(); ++i)
      {
        // modify sequence
        String mod_seq = res_seq;
        if (mod_seq[pos] != targets[i])
        {
          mod_seq[pos] = targets[i];
          generateTargetSequences(mod_seq, pos + 1, map_source2target, target_sequences);
        }
      }
      ++pos;
    }
  }

  // check and add only valid sequences (containing only target nucleotides or nucleotides that are both source and target nucleotides)
  Size count = 0;
  for (Size pos = 0; pos != res_seq.size(); ++pos)
  {
    TConstMapIterator target_iterator = map_source2target.find(res_seq[pos]);
    if (target_iterator == map_source2target.end()) // no pure source nucleotide?
    {
      count++;
    }
    else // check if source nucleotide is also a valid target nucleotide
    {
      const vector<char>& targets = target_iterator->second;
      for (Size i = 0; i != targets.size(); ++i)
      {
        if (res_seq[pos] == targets[i])
        {
          count++;
        }
      }
    }
  }

  if (count == res_seq.size())
  {
    target_sequences.push_back(res_seq);
  }
}
}

