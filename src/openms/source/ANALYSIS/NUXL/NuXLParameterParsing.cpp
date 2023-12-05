// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

// extract all possible marker ions and make them unique by mass
std::vector<NuXLFragmentAdductDefinition> NuXLParameterParsing::getMarkerIonsMassSet(const NuXLParameterParsing::PrecursorsToMS2Adducts& pc2adducts)
{  
  std::vector<NuXLFragmentAdductDefinition> marker_ions_unique_by_mass;
  for (const auto& p2a : pc2adducts)
  { // for all precursor adducts (e.g.: "UU-H2O") to all chemically feasible fragment adducts
    for (const NuXLFragmentAdductDefinition& f : p2a.second.marker_ions)
    { // add all marker ions
      marker_ions_unique_by_mass.push_back(f);
    }
  }
  // NuXLFragmentAdductDefinition has formula, name and mass. 
  // We are only interested in unique masses for spectra generation so we need to make the.

  auto less_by_mass = [](const NuXLFragmentAdductDefinition & a, const NuXLFragmentAdductDefinition & b) -> bool
    { 
      return a.mass < b.mass; 
    };

  sort(marker_ions_unique_by_mass.begin(), marker_ions_unique_by_mass.end(), less_by_mass);
  marker_ions_unique_by_mass.erase(
    unique(marker_ions_unique_by_mass.begin(), marker_ions_unique_by_mass.end(), less_by_mass), 
    marker_ions_unique_by_mass.end());

  return marker_ions_unique_by_mass;
}

NuXLParameterParsing::PrecursorsToMS2Adducts
NuXLParameterParsing::getAllFeasibleFragmentAdducts(
  const NuXLModificationMassesResult & precursor_adducts,
  const NuXLParameterParsing::NucleotideToFragmentAdductMap & nucleotide_to_fragment_adducts,
  const set<char> & can_xl, 
  const bool always_add_default_marker_ions,
  const bool default_marker_ions_RNA)
{
  PrecursorsToMS2Adducts all_pc_all_feasible_adducts;

  map<String, double> pc2mass;
  map<String, String> pc2ef;

  using PrecursorAdductMassAndMS2Fragments = pair<double, set<double>>;
  using PrecursorAdductAndXLNucleotide = pair<String, char>;
  map<PrecursorAdductMassAndMS2Fragments, vector<PrecursorAdductAndXLNucleotide>> mass_frags2pc_xlnuc;
  
  // for all distinct precursor adduct formulas/masses
  for (auto const & pa : precursor_adducts.formula2mass)
  {
    // get all precursor nucleotide formulas matching current empirical formula/mass
    const String& ef = pa.first;
    const double pc_mass = pa.second;
    const auto& ambiguities = precursor_adducts.mod_combinations.at(ef);

    if (ambiguities.size() >= 2)
    {
      OPENMS_LOG_DEBUG << ambiguities.size() << " nucleotide formulas are ambiguous on the level of empirical formula: " << endl;
      for (auto const & pc_adduct : ambiguities) { OPENMS_LOG_DEBUG << pc_adduct << endl; }
    } 

    // for each ambiguous (at the level of empirical formula) precursor adduct (stored as nucleotide formula e.g.: "AU-H2O")
    for (const String & pc_adduct : ambiguities)
    {
      // calculate feasible fragment adducts and store them for lookup
      const MS2AdductsOfSinglePrecursorAdduct& feasible_adducts = getFeasibleFragmentAdducts(pc_adduct, 
        ef, 
        nucleotide_to_fragment_adducts, 
        can_xl,
        always_add_default_marker_ions,
        default_marker_ions_RNA);

      // TODO: check if needed anymore - std::sort(feasible_adducts.begin(), feasible_adducts.end());
       all_pc_all_feasible_adducts[pc_adduct] = feasible_adducts;
       pc2mass[pc_adduct] = pc_mass;
       pc2ef[pc_adduct] = ef;
    }
  }

  // print feasible fragment adducts and marker ions
  for (auto const & fa : all_pc_all_feasible_adducts)
  {
    const String & pc = fa.first;
    OPENMS_LOG_DEBUG << "Precursor adduct: " << pc << "\t" << pc2ef[pc] << "\t" << pc2mass[pc] << "\n";

    // collect set of masses to detect ambiguities
    set<double> mi_fragments;
    for (auto const & mis : fa.second.marker_ions) { mi_fragments.insert(mis.mass); }

    for (auto const & ffa : fa.second.feasible_adducts)
    {
      const char & nucleotide = ffa.first;
      set<double> fragments;
      for (auto const & f : mi_fragments) { fragments.insert(f); } // copy over common marker ion masses
      for (auto const & f : ffa.second) { fragments.insert(f.mass); } // copy over xl-ed nucleotide associated fragment masses
      // we want to track if a certain MS1 mass and set of MS2 ions allow to distinguish 
      // (at least in theory) precursor adduct and cross-linked nucleotide 
      // to do so, we map mass and set of fragment ions to precursor adduct and cross-linked nucleotide 
      mass_frags2pc_xlnuc[{pc2mass[pc], fragments}].push_back({pc, nucleotide});

      OPENMS_LOG_DEBUG << "  Cross-linkable nucleotide '" << nucleotide << "' and feasible fragment adducts:" << endl;
      for (auto const & a : ffa.second)
      {
        OPENMS_LOG_DEBUG << "    " << a.name << "\t" << a.formula.toString() << "\t" << a.mass << "\n";
      }
    }

    OPENMS_LOG_DEBUG << "  Marker ions." << endl;
    for (auto const & ffa : fa.second.marker_ions)
    {
      OPENMS_LOG_DEBUG << "    "  << ffa.name << "\t" << ffa.formula.toString() << "\t" << ffa.mass << "\n";
    }
  }
  OPENMS_LOG_DEBUG << endl;

  // report precursor adducts/cross-linked nucleotide combinations with same mass that produce the exact same MS2 ions as indistinguishable
  // e.g., for nucleotides that only produce ribose fragments there might be a lot of overlap
  for (auto const & f : mass_frags2pc_xlnuc)
  {
    if (f.second.size() > 1)
    {
      OPENMS_LOG_DEBUG << "Theoretical indistinguishable cross-link adducts detected: "  << endl;
      for (auto const & m : f.second)
      {
        OPENMS_LOG_DEBUG << "\tPrecursor: " << m.first << "\t and cross-linked nucleotide: " << m.second << endl;
        for (auto const & fragment_mass : f.first.second)
        {
          OPENMS_LOG_DEBUG << "\tFragment mass: " << fragment_mass << endl;
        }
      }      
    }
  }

  return all_pc_all_feasible_adducts;
}

NuXLParameterParsing::NucleotideToFragmentAdductMap
NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(StringList fragment_adducts)
{
  NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts;
  for (String t : fragment_adducts)
  {
    t.removeWhitespaces();

    EmpiricalFormula formula;
    String name;

    // format is: target_nucletide:formula;name
    char target_nucleotide = t[0];
    if (t[1] != ':')
    {
      OPENMS_LOG_WARN << "Missing ':'. Wrong format of fragment_adduct string: " << t << endl;
      return NucleotideToFragmentAdductMap();
    }

    // remove target nucleotide and : character from t
    t = t.substr(2);

    // split into formula and name
    vector<String> fs;
    t.split(";", fs);
    if (fs.size() == 1) // no name provided so we just take the formula as name
    {
      formula = EmpiricalFormula(fs[0]);
      name = formula.toString();
    }
    else if (fs.size() == 2)
    {
      formula = EmpiricalFormula(fs[0]);
      name = fs[1];
    }
    else
    {
      OPENMS_LOG_WARN << "Wrong format of fragment_adduct string: " << t << endl;
      return NucleotideToFragmentAdductMap();
    }

    NuXLFragmentAdductDefinition fad;
    fad.name = name;
    fad.formula = formula;
    fad.mass = formula.getMonoWeight();

    nucleotide_to_fragment_adducts[target_nucleotide].insert(fad);

    // register all fragment adducts as N- and C-terminal modification (if not already registered)
    if (!ModificationsDB::getInstance()->has(name))
    {
      std::unique_ptr<ResidueModification> c_term{new ResidueModification()};
      c_term->setId(name);
      c_term->setName(name);
      c_term->setFullId(name + " (C-term)");
      c_term->setTermSpecificity(ResidueModification::C_TERM);
      c_term->setDiffMonoMass(fad.mass);
      ModificationsDB::getInstance()->addModification(std::move(c_term));

      std::unique_ptr<ResidueModification> n_term{new ResidueModification()};
      n_term->setId(name);
      n_term->setName(name);
      n_term->setFullId(name + " (N-term)");
      n_term->setTermSpecificity(ResidueModification::N_TERM);
      n_term->setDiffMonoMass(fad.mass);
      ModificationsDB::getInstance()->addModification(std::move(n_term));
    }
  }

#ifdef DEBUG_OpenNuXL
  for (auto const & p2fas : precursor_to_fragment_adducts)
    {
      for (auto const & p2fa : p2fas.second)
      {
        cout << "nucleotide:" << p2fas.first
             << " fragment adduct:" << p2fa.formula.toString()
             << " fragment adduct mass:" << p2fa.mass
             << " name:" <<  p2fa.name << endl;
      }
    }
#endif

  return nucleotide_to_fragment_adducts;
}




MS2AdductsOfSinglePrecursorAdduct
NuXLParameterParsing::getFeasibleFragmentAdducts(const String &exp_pc_adduct,
                                                  const String &exp_pc_formula,
                                                  const NuXLParameterParsing::NucleotideToFragmentAdductMap &nucleotide_to_fragment_adducts,
                                                  const set<char> &can_xl,
                                                  const bool always_add_default_marker_ions,
                                                  const bool default_marker_ions_RNA)
{
  OPENMS_LOG_DEBUG << "Generating fragment adducts for precursor adduct: '" << exp_pc_adduct << "'" << endl;

  MS2AdductsOfSinglePrecursorAdduct ret;

  // no precursor adduct? 
  if (exp_pc_formula.empty()) { return ret; } // no fragment adducts or marker ions are expected!

  // count nucleotides in precursor adduct (e.g.: "TCA-H2O" yields map: T->1, C->1, A->1)
  // and determine the set of cross-linkable nucleotides in the precursor adduct
  size_t nt_count(0);

  set<char> exp_pc_xl_nts; // the cross-linkable nucleotides in the precursor adduct
  map<char, Size> exp_pc_nucleotide_count; // all nucleotides in the precursor adduct (e.g., used to determine marker ions)

  {
    for (String::const_iterator exp_pc_it = exp_pc_adduct.begin(); exp_pc_it != exp_pc_adduct.end(); ++exp_pc_it, ++nt_count)
    {
      // we are finished with nucleotides in string if first loss/gain is encountered
      if (*exp_pc_it == '+' || *exp_pc_it == '-') break;

      // count occurence of nucleotide
      if (exp_pc_nucleotide_count.count(*exp_pc_it) == 0)
      {
        exp_pc_nucleotide_count[*exp_pc_it] = 1;
        if (can_xl.count(*exp_pc_it)) { exp_pc_xl_nts.insert(*exp_pc_it); };
      }
      else
      {
        exp_pc_nucleotide_count[*exp_pc_it]++;
      }
    }
  }

  // check if at least one nucleotide present that can cross link
  bool has_xl_nt = !exp_pc_xl_nts.empty();
  OPENMS_LOG_DEBUG << "\t" << exp_pc_adduct << " has cross-linkable nucleotide (0 = false, 1 = true): " << has_xl_nt << endl;

  // no cross-linkable nt contained in the precursor adduct? Return an empty fragment adduct definition set
  if (!has_xl_nt) { return ret; }

  // determine if there is a nucleotide/sugar/etc. that must be the cross-linked one
  set<char> must_xl;
  for (auto c : exp_pc_xl_nts) { if (c == 'd' || c == 'r') must_xl.insert(c); }
  if (must_xl.size() >= 2) 
  {
    OPENMS_LOG_WARN << "More than one nucleotide present that is marked as mandatory cross-linked ('d' or 'r')." << endl;
    return ret; 
  }
  else if (must_xl.size() == 1)
  {  
    cout << "Mandatory cross-linking nt/sugar: " << *must_xl.begin() << " in precursor adduct: " << exp_pc_adduct  << endl;
    exp_pc_xl_nts = must_xl;
  } 
  // else we have no mandatory cross-linked nts

  ///////////////////////////////////////////////////////////////////
  // HERE: at least one cross-linkable nt present in precursor adduct

  // extract loss string from precursor adduct (e.g.: "-H2O")
  // String exp_pc_loss_string(exp_pc_it, exp_pc_adduct.end());

  OPENMS_LOG_DEBUG << "\t" << exp_pc_adduct << " is monomer (1 = true, >1 = false): " << nt_count << endl;

  // Handle the cases of monomer or oligo nucleotide bound to the precursor.
  // This distinction is made because potential losses on the precursor only allows us to reduce the set of chemical feasible fragment adducts if they are on a monomer.
  // In the case of an oligo we can't be sure if the cross-linked amino acid or any other in the oligo had the loss.
  if (nt_count > 1)  // No monomer? For every nucleotide that can be cross-linked: Create all fragment adducts without restriction by losses (no restriction as the loss could be on the other nts)
  {
    // for each nucleotide and potential set of fragment adducts
    for (auto const & n2fa : nucleotide_to_fragment_adducts)
    {
      const char & nucleotide = n2fa.first; // the nucleotide without any associated loss
      const set<NuXLFragmentAdductDefinition>& fragment_adducts = n2fa.second; // all potential fragment adducts that may arise from the unmodified nucleotide

      // check if nucleotide is cross-linkable and part of the precursor adduct
      if (exp_pc_xl_nts.find(nucleotide) != exp_pc_xl_nts.end())
      {
        OPENMS_LOG_DEBUG << "\t" << exp_pc_adduct << " found nucleotide: " << String(nucleotide) << " in precursor RNA." << endl;
        OPENMS_LOG_DEBUG << "\t" << exp_pc_adduct << " nucleotide: " << String(nucleotide) << " has fragment_adducts: " << fragment_adducts.size() << endl;

        // store feasible adducts associated with a cross-link with character nucleotide
        vector<NuXLFragmentAdductDefinition> faa;
        std::copy(fragment_adducts.begin(), fragment_adducts.end(), back_inserter(faa));
        ret.feasible_adducts.emplace_back(make_pair(nucleotide, faa));
      }
    }
    
    // Create set of marker ions for *all* nucleotides contained in the precursor adduct (including those that do not cross-link.)
    // Note: The non-cross-linked nt in the precursor adduct are more likely to produce the marker ions (=more fragile than the cross-linked nt).
    set<NuXLFragmentAdductDefinition> marker_ion_set;
    for (auto const & n2fa : nucleotide_to_fragment_adducts)
    {
      const char & nucleotide = n2fa.first; // the nucleotide without any associated loss
      if (exp_pc_nucleotide_count.find(nucleotide) != exp_pc_nucleotide_count.end())
      { 
        marker_ion_set.insert(n2fa.second.begin(), n2fa.second.end()); 
      }
    }
    std::move(std::begin(marker_ion_set), std::end(marker_ion_set), std::back_inserter(ret.marker_ions));
  }
  else // nt_count == 1: monomer. We need to check if the neutral loss reduces the set of feasible (e.g., chemically sound) fragment adducts
  {
    for (auto const & n2fa : nucleotide_to_fragment_adducts)
    {
      const char & nucleotide = n2fa.first; // one letter code of the nt

      set<NuXLFragmentAdductDefinition> fas = n2fa.second; // all potential fragment adducts that may arise from NT assuming no loss on the precursor adduct

      // check if nucleotide is cross-linkable and part of the precursor adduct
      if (exp_pc_xl_nts.find(nucleotide) != exp_pc_xl_nts.end())
      {
        OPENMS_LOG_DEBUG << "\t" << exp_pc_adduct << " found nucleotide: " << String(nucleotide) << " in precursor NA adduct." << endl;
        OPENMS_LOG_DEBUG << "\t" << exp_pc_adduct << " nucleotide: " << String(nucleotide) << " has fragment_adducts: " << fas.size() << endl;

        // check chemical feasibility by checking if subtraction of adduct would result in negative elemental composition
        for (auto it = fas.begin(); it != fas.end(); )
        {
          bool negative_elements = (EmpiricalFormula(exp_pc_formula) - it->formula).toString().hasSubstring("-");

          if (negative_elements) // fragment adduct can't be subformula of precursor adduct
          {
            it = fas.erase(it);
          }
          else
          {
            ++it; // STL erase idiom (mind the pre-increment)
          }
        }
        
        if (fas.empty()) continue; // no feasible fragment adducts left? continue

        // store feasible adducts associated with a cross-link with character nucleotide[0]
        vector<NuXLFragmentAdductDefinition> faa;
        std::copy(fas.begin(), fas.end(), back_inserter(faa));
        ret.feasible_adducts.emplace_back(make_pair(nucleotide, faa));

        // We only have one nucleotide in the precursor adduct (the cross-linked one)
        // Note: marker ions of the cross-linked nucleotide are often missing or of very low intensity
        std::copy(std::begin(fas), std::end(fas), std::back_inserter(ret.marker_ions));
      }
    }
  }

  // for chemical cross-linkers like DEB, fragments always carry DEB but marker ions might just be U or U' (and losses)
  // In that case we need to ensure that the the default marker ions are added
  if (always_add_default_marker_ions)
  {
    // Note: add the uncharged mass. Protons are added during spectrum generation.
    if (default_marker_ions_RNA) // TODO: check if we can derive this from target nucleotides
    {
      ret.marker_ions.push_back({EmpiricalFormula("C9H13N2O9P1"), String("U"), EmpiricalFormula("C9H13N2O9P1").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C9H14N3O8P"), "C", EmpiricalFormula("C9H14N3O8P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C10H14N5O8P"), "G", EmpiricalFormula("C10H14N5O8P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C10H14N5O7P"), "A", EmpiricalFormula("C10H14N5O7P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C4H4N2O2"), "U'", EmpiricalFormula("C4H4N2O2").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C4H5N3O"), "C'", EmpiricalFormula("C4H5N3O").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C5H5N5O"), "G'",  EmpiricalFormula("C5H5N5O").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C5H5N5"), "A'", EmpiricalFormula("C5H5N5").getMonoWeight()});
    }
    else // DNA
    {
      ret.marker_ions.push_back({EmpiricalFormula("C10H15N2O8P"), "T", EmpiricalFormula("C10H15N2O8P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C9H14N3O7P"), "C", EmpiricalFormula("C9H14N3O7P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C10H14N5O7P"), "G", EmpiricalFormula("C10H14N5O7P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C10H14N5O6P"), "A", EmpiricalFormula("C10H14N5O6P").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C5H6N2O2"), "T'", EmpiricalFormula("C5H6N2O2").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C4H5N3O"), "C'", EmpiricalFormula("C4H5N3O").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C5H5N5O"), "G'", EmpiricalFormula("C5H5N5O").getMonoWeight()});
      ret.marker_ions.push_back({EmpiricalFormula("C5H5N5"), "A'", EmpiricalFormula("C5H5N5").getMonoWeight()});
    } 
  }


  // Because, e.g., ribose might be a feasible fragment of any nucleotide, we keep only one version
  // Note: sort by formula and (as tie breaker) the name
  std::sort(ret.marker_ions.begin(), ret.marker_ions.end(),
    [](NuXLFragmentAdductDefinition const & a, NuXLFragmentAdductDefinition const & b)
    {
      const String as = a.formula.toString();
      const String bs = b.formula.toString();
      return std::tie(as, a.name) < std::tie(bs, b.name);
    }
  );
  // Note: for uniqueness, we only rely on the formula (in case of tie: keeping the first = shortest name)
  auto it = std::unique(ret.marker_ions.begin(), ret.marker_ions.end(),
    [](NuXLFragmentAdductDefinition const & a, NuXLFragmentAdductDefinition const & b)
    {
      return a.formula == b.formula;
    }
  );
  ret.marker_ions.resize(std::distance(ret.marker_ions.begin(), it));

  // print feasible fragment adducts
  for (auto const & ffa : ret.feasible_adducts)
  {
    const char & nucleotide = ffa.first;
    OPENMS_LOG_DEBUG << "  Cross-linkable nucleotide '" << nucleotide << "' and feasible fragment adducts:" << endl;
    for (auto const & a : ffa.second)
    {
      OPENMS_LOG_DEBUG << "\t" << a.name << "\t" << a.formula.toString() << "\t" << a.mass << "\n";
    }
  }

  // print marker ions
  OPENMS_LOG_DEBUG << "  Marker ions:" << endl;
  for (auto const & a : ret.marker_ions)
  {
    OPENMS_LOG_DEBUG << "\t" << a.name << "\t" << a.formula.toString() << "\t" << a.mass << "\n";
  }

  return ret;
}


vector<ResidueModification> NuXLParameterParsing::getModifications(StringList modNames) 
{
  vector<ResidueModification> modifications;

  // iterate over modification names and add to vector
  for (String modification : modNames)
  {
    ResidueModification rm;
    if (modification.hasSubstring(" (N-term)"))
    {
      modification.substitute(" (N-term)", "");
      rm = *ModificationsDB::getInstance()->getModification(modification, "", ResidueModification::N_TERM);
    }
    else if (modification.hasSubstring(" (C-term)"))
    {
      modification.substitute(" (C-term)", "");
      rm = *ModificationsDB::getInstance()->getModification(modification, "", ResidueModification::C_TERM);
    }
    else
    {
      rm = *ModificationsDB::getInstance()->getModification(modification);
    }
    modifications.push_back(rm);
  }

  return modifications;
}


}

