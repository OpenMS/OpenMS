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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{

//static
bool RNPxlModificationsGenerator::notInSeq(String res_seq, String query)
{
  // special case: empty query is in every seq -> false
  if (query.empty()) { return false; }

  // test all k-mers with k=size of query
  for (Int l = 0; l <= (Int)res_seq.size() - (Int)query.size(); ++l)
  {
    String a = res_seq.substr(l, query.size());
    String b = query;

    sort(a.begin(), a.end());
    sort(b.begin(), b.end());

    if (a == b) { return false; }
  }
  return true;
}

//static
RNPxlModificationMassesResult RNPxlModificationsGenerator::initModificationMassesRNA(StringList target_nucleotides,
                                                                                     StringList nt_groups,
                                                                                     std::set<char> can_xl,
                                                                                     StringList mappings,
                                                                                     StringList modifications,
                                                                                     String sequence_restriction,
                                                                                     bool cysteine_adduct,
                                                                                     Int max_length)
{
  String original_sequence_restriction = sequence_restriction;

  // 152 modification
  const String cysteine_adduct_string("C4H8S2O2");//FIXME: why is this changed from ancestor?
  const EmpiricalFormula cysteine_adduct_formula(cysteine_adduct_string); // 152 modification

  RNPxlModificationMassesResult result;

  // read nucleotides and empirical formula of monophosphate
  map<String, EmpiricalFormula> map_target_to_formula;
  for (auto const & s : target_nucleotides)
  {
    vector<String> fields;
    s.split("=", fields);
    map_target_to_formula[fields[0]] = EmpiricalFormula(fields[1]);
  }

  // read mapping of source to target
  map<char, vector<char> > map_source_to_targets;
  for (auto const & s : mappings)
  {
    vector<String> fields;
    s.split("->", fields);
    map_source_to_targets[fields[0][0]].push_back(fields[1][0]);
  }

  // extract source nucleotides based on mapping
  vector<char> source_nucleotides; // nucleotides as expected in the restriction sequence
  for (auto const & s : mappings)
  {
    source_nucleotides.push_back(s[0]);
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

  // erase trivial cases from mapping so only the combinatorial cases: 1 source -> n targets remain
  for (map<char, vector<char> >::iterator sit = map_source_to_targets.begin(); sit != map_source_to_targets.end(); )
  {
    char source = sit->first;
    char first_target = sit->second[0];

    if (sit->second.size() == 1 && source == first_target) // trivial case e.g. A->A... no substitution needed
    {
      map_source_to_targets.erase(sit++);
    }
    else if (sit->second.size() == 1 && source != first_target) // simple rename e.g. A->X... simply substitute all in restriction sequence
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
    OPENMS_LOG_WARN << "WARNING: no restriction on sequence but multiple target nucleotides specified. Will generate huge amount of sequences" << endl;
  }

  using NucleotideModificationSubFormula = pair<EmpiricalFormula, bool>; // e.g., "H2O", true
  using NucleotideModification = vector<NucleotideModificationSubFormula>;
  using NucleotideModifications = vector<NucleotideModification>;

  // map nucleotide to list of empirical MS1 precursor losses/gains
  // nucleotide->all loss/gain formulas (each composed of subformulas)->subformulas
  map<String, NucleotideModifications> map_to_nucleotide_modifications;

  for (String m : modifications)
  {
    // extract target nucleotide
    if (m[1] != ':') 
    { 
      throw Exception::MissingInformation(
       __FILE__, 
       __LINE__, 
       OPENMS_PRETTY_FUNCTION, 
       " Modifications parameter must specify nucleotide and forulas in format 'U:+H2O-H2O'.");
    };

    String target_nucleotide = m[0];

    NucleotideModification nucleotide_modification;

    m = m.substr(2); // remove nucleotide and ':' from front of string

    // decompose string into subformulas
    m.substitute("-", "#-");
    m.substitute("+", "#+");
    vector<String> ems;
    m.split("#", ems);
    for (Size j = 0; j != ems.size(); ++j)
    {
      if (ems[j].empty()) { continue; }

      bool mod_is_subtractive(false);

      if (ems[j][0] == '-')
      {
        mod_is_subtractive = true;
        ems[j].remove('-');
      }
      else if (ems[j][0] == '+')
      {
        ems[j].remove('+');
      }

      EmpiricalFormula ef(ems[j]);
      ef.setCharge(0);
      NucleotideModificationSubFormula sf = make_pair(ef, mod_is_subtractive);
      nucleotide_modification.push_back(sf);
    }
    // add formula to target nucleotide
    map_to_nucleotide_modifications[target_nucleotide].push_back(nucleotide_modification);
  }

  // generate all target sequences by substituting each source nucleotide by their target nucleotide(s)
  StringList target_sequences;
  generateTargetSequences(sequence_restriction, 0, map_source_to_targets, target_sequences);

  OPENMS_LOG_INFO << "target sequence(s):" << target_sequences.size() << endl;

  if (!original_sequence_restriction.empty())
  {
    for (Size i = 0; i != target_sequences.size(); ++i)
    {
      if (target_sequences[i].size() < 60)
      {
        OPENMS_LOG_INFO << target_sequences[i] << endl;
      }
      else
      {
        OPENMS_LOG_INFO << target_sequences[i].prefix(60) << "..."  << endl;
      }
    }
  }

  {
    // Append precursor modifications (e.g., "-H2O") 
    // to generate modified nucleotides: e.g.: "U" -> "U", "U-H2O", ...
    vector<EmpiricalFormula> actual_combinations;
    for (map<String, EmpiricalFormula>::const_iterator mit = map_target_to_formula.begin(); mit != map_target_to_formula.end(); ++mit)
    {
      String target_nucleotide = mit->first;
      OPENMS_LOG_INFO << "target nucleotide: " << target_nucleotide << endl;

      EmpiricalFormula target_nucleotide_formula = mit->second;

      NucleotideModifications nt_mods = map_to_nucleotide_modifications[target_nucleotide];

      for (NucleotideModification const & nt_mod : nt_mods) // loop over list of nucleotide specific modifications
      {
        EmpiricalFormula e(target_nucleotide_formula);
        String s(target_nucleotide);
        for (NucleotideModificationSubFormula const & sf : nt_mod) // loop over subformulae
        {
          // concatenate additive / subtractive substrings (e.g., "+H2O", "-H3PO")
          EmpiricalFormula mod_ef(sf.first);
          String mod(sf.first.toString());
          if (sf.second) 
          {  // subtractive
            e = e - mod_ef;
           s += "-" + mod;
           }
          else 
          {  // additive
            e = e + mod_ef;
            s += "+" + mod;;
          }
        }
        actual_combinations.push_back(e);
        result.mod_combinations[actual_combinations.back().toString()].insert(s);
        OPENMS_LOG_INFO << "\t" << "modifications: " << s << "\t\t" << e.toString() << endl;
      }
    }

    // Generate >=1 nucleotide precursor adducts (e.g., "UU-H2O-H3PO")
    // In every loop iteration, an unmodified target_nucleotide (e.g., "U", "A", ... ) is added to the chain
    // The first element of the chain is an unmodified AND modified nucleotides.
    // That way, at most one modified nucleotide is part of the chain
    vector<EmpiricalFormula> all_combinations = actual_combinations;
    for (Int i = 0; i < max_length - 1; ++i)
    {
      vector<EmpiricalFormula> new_combinations;
      for (auto const & target_to_formula : map_target_to_formula) // loop over nucleotides (unmodified)
      {
        const String & target_nucleotide = target_to_formula.first;
        const EmpiricalFormula & target_nucleotide_formula = target_to_formula.second;

        for (EmpiricalFormula const & ac : actual_combinations) // append unmodified nucleotide to yield a (i+1)-mer
        {
          new_combinations.push_back(target_nucleotide_formula + ac - EmpiricalFormula("H2O")); // -H2O because of condensation reaction
          all_combinations.push_back(target_nucleotide_formula + ac - EmpiricalFormula("H2O")); // " "
          const set<String>& ambiguities = result.mod_combinations[ac.toString()];
          for (auto const & s : ambiguities)
          {
            result.mod_combinations[all_combinations.back().toString()].insert(target_nucleotide + s);
            OPENMS_LOG_DEBUG << target_nucleotide + s << endl;
          }
        }
      }
      actual_combinations = new_combinations;
    }

    for (Size i = 0; i != all_combinations.size(); ++i)
    {
      result.mod_masses[all_combinations[i].toString()] = all_combinations[i].getMonoWeight();
    }
  }

  std::cout << "Filtering on restrictions... " << endl;

  // Remove precursor adducts that
  // 1) do not contain a cross-linkable nucleotide
  // 2) or contain no cross-linkable nucleotide that is part of the restricted target sequences
  std::vector<pair<String, String> > violates_restriction; // elemental composition, nucleotide style formula
  for (Map<String, double>::ConstIterator mit = result.mod_masses.begin(); mit != result.mod_masses.end(); ++mit)
  {
    // remove additive or subtractive modifications from string as these are not used in string comparison
    const set<String>& ambiguities = result.mod_combinations[mit->first];
    for (String const & s : ambiguities)
    {
      String nucleotide_style_formula(s);

      // get nucleotide formula without losses / gains (e.g., "U" instead of "U-H2O")
      Size p1 = nucleotide_style_formula.find('-');
      Size p2 = nucleotide_style_formula.find('+');
      Size p = min(p1, p2);
      if (p != String::npos)
      {
        nucleotide_style_formula = nucleotide_style_formula.prefix(p);
      }

      // check if nucleotide formula contains a cross-linkable amino acid
      bool has_xl_nt(false);
      for (auto const & c : nucleotide_style_formula) { if (can_xl.count(c) > 0) { has_xl_nt = true; break;};  }

      if (!has_xl_nt) 
      { // no cross-linked nucleotide => not valid
        violates_restriction.push_back(make_pair(mit->first, s)); 
        continue;
      }

      // check if nucleotides from more than one nt_group are present (e.g. from DNA and RNA)
      Size found_in_n_groups(0);
      for (const String & n : nt_groups)
      { 
        if (nucleotide_style_formula.find_first_of(n) != string::npos) { ++found_in_n_groups; }
      }
      // nucleotide stile formula (e.g. AATU matches to more than one group (e.g., RNA and DNA))?
      if (found_in_n_groups > 1)
      {
        violates_restriction.push_back(make_pair(mit->first, s)); 
        continue;
      }

      // check if nucleotide is contained in at least one of the target sequences
      bool containment_violated(false);
      Size violation_count(0);
      for (auto const & current_target_seq : target_sequences)
      {
        if (notInSeq(current_target_seq, nucleotide_style_formula)) { ++violation_count; }
      }

      if (violation_count == target_sequences.size()) { containment_violated = true; }

      if (containment_violated)
      { 
        violates_restriction.push_back(make_pair(mit->first, s)); // chemical formula, nucleotide style formula pair violates restrictions
      }
    }
  }

  for (size_t i = 0; i != violates_restriction.size(); ++i)
  {
    const String& chemical_formula = violates_restriction[i].first;
    result.mod_combinations[chemical_formula].erase(violates_restriction[i].second);
    OPENMS_LOG_DEBUG << "filtered sequence: " 
              << chemical_formula 
              << "\t" 
              << violates_restriction[i].second << endl;
  }

  // standard associative-container erase idiom
  for (map<String, set<String> >::iterator mcit = result.mod_combinations.begin(); mcit != result.mod_combinations.end(); )
  {
    if (mcit->second.empty())
    {
      result.mod_masses.erase(mcit->first); // remove from mod masses
      result.mod_combinations.erase(mcit++); // don't change precedence !
    }
    else
    {
      ++mcit;   // don't change precedence !
    }
  }


  // Optional: add cystein (DTT) adduct
  if (cysteine_adduct)
  {
    result.mod_masses[cysteine_adduct_formula.toString()] = cysteine_adduct_formula.getMonoWeight();
    result.mod_combinations[cysteine_adduct_formula.toString()].insert(cysteine_adduct_string);
  }

  // output index -> empirical formula -> (ambiguous) nucleotide formulas
  // nucleotide formulas which only differ in nucleotide ordering are only printed once
  // e.g. 5 C19H24N7O12P1 573.122 ( AU-H1O3P1 )
  double pseudo_rt = 1;
  for (auto const & m : result.mod_masses)
  {
    if (cysteine_adduct && m.first == cysteine_adduct_formula.toString())
    {
      OPENMS_LOG_INFO << "Precursor adduct " << pseudo_rt++ << "\t:\t" << m.first << " " << m.second << " ( cysteine adduct )" << endl;
      continue;
    }

    OPENMS_LOG_INFO << "Precursor adduct " << pseudo_rt++ << "\t:\t" << m.first << " " << m.second << " ( ";

    const set<String>& ambiguities = result.mod_combinations[m.first];
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


      // only print ambiguous sequences once
      if (printed.find(nucleotide_style_formula) == printed.end())
      {
        OPENMS_LOG_INFO << nucleotide_style_formula << " ";
        printed.insert(nucleotide_style_formula);
      }
      else
      {
        OPENMS_LOG_DEBUG << "ambigious " << nucleotide_style_formula << endl;
      }
    }

    OPENMS_LOG_INFO << ")" << endl;
  }
  OPENMS_LOG_INFO << "Finished generation of modification masses." << endl;
  return result;
}

//static
void  RNPxlModificationsGenerator::generateTargetSequences(const String& res_seq,
                                                           Size param_pos,
                                                           const map<char, vector<char> >& map_source2target,
                                                           StringList& target_sequences)
{
  while (param_pos < res_seq.size())
  {
    // check if current character is in source 2 target map
    auto target_iterator = map_source2target.find(res_seq[param_pos]);
    if (target_iterator == map_source2target.end())
    {
      ++param_pos;
    }
    else // yes?
    {
      const vector<char>& targets = target_iterator->second;
      for (Size i = 0; i != targets.size(); ++i)
      {
        // modify sequence
        String mod_seq = res_seq;
        if (mod_seq[param_pos] != targets[i])
        {
          mod_seq[param_pos] = targets[i];
          generateTargetSequences(mod_seq, param_pos + 1, map_source2target, target_sequences);
        }
      }
      ++param_pos;
    }
  }

  // check and add only valid sequences (containing only target nucleotides or nucleotides that are both source and target nucleotides)
  Size count = 0;
  for (Size pos = 0; pos != res_seq.size(); ++pos)
  {
    auto target_iterator = map_source2target.find(res_seq[pos]);

    // no pure source nucleotide?
    if (target_iterator == map_source2target.end())
    {
      count++;
    }
    else // check if source nucleotide is also a valid target nucleotide
    {
      const vector<char>& targets = target_iterator->second;
      for (Size i = 0; i != targets.size(); ++i)
      {
        if (res_seq[pos] == targets[i]) { count++; }
      }
    }
  }

  if (count == res_seq.size())
  {
    target_sequences.push_back(res_seq);
  }
}

}

