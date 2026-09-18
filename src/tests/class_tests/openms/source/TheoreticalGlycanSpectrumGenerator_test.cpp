// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/MzPAF.h>
#include <OpenMS/CHEMISTRY/TheoreticalGlycanSpectrumGenerator.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/test_config.h>
#include <algorithm>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tuple>

using namespace OpenMS;
using Generator = TheoreticalGlycanSpectrumGenerator;
using Composition = Generator::Composition;
using Ion = Generator::IonType;

namespace
{
Composition comp(std::initializer_list<std::pair<std::string, int>> counts)
{
  Composition result;
  for (const auto& [symbol, count] : counts)
  {
    result.components.emplace_back(symbol, count);
  }
  return result;
}

const Generator::Fragment& findIon(const std::vector<Generator::Fragment>& fragments, std::string_view name, Int charge = 1)
{
  const auto it = std::find_if(fragments.begin(), fragments.end(), [&](const auto& f) { return f.name == name && f.charge == charge; });
  if (it == fragments.end()) { throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Missing ion: " + std::string(name)); }
  return *it;
}

using GlypyKey = std::tuple<Ion, std::optional<Size>, std::vector<Size>>;
struct GlypyCase
{
  GlycanStructure tree;
  std::map<GlypyKey, double> masses;
};

std::map<std::string, GlypyCase> loadGlypyCases()
{
  std::ifstream input(OPENMS_GET_TEST_DATA_PATH("TheoreticalGlycanSpectrumGenerator_glypy/structural.tsv"));
  if (! input) { throw std::runtime_error("Cannot open GlyPy structural fixture"); }
  std::map<std::string, GlypyCase> cases;
  const std::map<char, Ion> ion_types {{'B', Ion::B}, {'C', Ion::C}, {'Y', Ion::Y}, {'Z', Ion::Z}};
  for (std::string line; std::getline(input, line);)
  {
    if (line.empty() || line[0] == '#') { continue; }
    std::istringstream row(line);
    char record;
    std::string fixture;
    row >> record >> fixture;
    auto& reference = cases[fixture];
    if (record == 'N')
    {
      Size index;
      Int parent;
      std::string symbol;
      row >> index >> parent >> symbol;
      if (! row) { throw std::runtime_error("Malformed GlyPy node: " + line); }
      const auto attachment = parent < 0 ? std::nullopt : std::optional<Size>(parent);
      const auto actual_index = symbol.starts_with("Formula:")
                                  ? reference.tree.addMonosaccharide(ProForma::FormulaTag {symbol.substr(8), std::nullopt}, attachment)
                                  : reference.tree.addMonosaccharide(symbol, attachment);
      if (actual_index != index) { throw std::runtime_error("Unexpected GlyPy node order: " + line); }
    }
    else if (record == 'F')
    {
      char series;
      Int root;
      std::string branches;
      double mass;
      row >> series >> root >> branches >> mass;
      if (! row) { throw std::runtime_error("Malformed GlyPy fragment: " + line); }
      std::vector<Size> cuts;
      if (branches != "-")
      {
        std::replace(branches.begin(), branches.end(), ',', ' ');
        std::istringstream branch_stream(branches);
        for (Size cut; branch_stream >> cut;)
        {
          cuts.push_back(cut);
        }
      }
      const auto root_cut = root < 0 ? std::nullopt : std::optional<Size>(root);
      const auto inserted = reference.masses.emplace(GlypyKey {ion_types.at(series), root_cut, cuts}, mass);
      if (! inserted.second) { throw std::runtime_error("Duplicate GlyPy interpretation: " + line); }
    }
    else { throw std::runtime_error("Unknown GlyPy record: " + line); }
  }
  return cases;
}
} // namespace

START_TEST(TheoreticalGlycanSpectrumGenerator, "$Id$")
TOLERANCE_ABSOLUTE(0.00001)
TOLERANCE_RELATIVE(1.00000001)

START_SECTION((diagnostic ions from a ProForma composition))
{
  Generator::Options options;
  options.add_b_ions = options.add_y_ions = false;
  Generator generator(options);
  const auto parsed = ProForma::parse("N[Glycan:HexNAc2Hex5NeuAc1]ST");
  const auto& residue = std::get<ProForma::SequenceElement>(parsed.sequence[0]);
  const auto glycan = std::get<Composition>(residue.modifications[0].alternatives[0].first);
  const auto ions = generator.getFragments(glycan);
  TEST_EQUAL(ions.size(), 12)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc").getMZ(), 204.08665)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-H2O").getMZ(), 186.07608)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-H4O2").getMZ(), 168.06552)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-C2H4O2").getMZ(), 144.06552)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-CH6O3").getMZ(), 138.05495)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc-C2H6O3").getMZ(), 126.05495)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:Hex1HexNAc1").getMZ(), 366.13947)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:Neu5Ac").getMZ(), 292.10269)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:Neu5Ac-H2O").getMZ(), 274.09213)
  const auto aliases = generator.getFragments(comp({{"NeuAc", 1}, {"Neu5Ac", 2}, {"NeuGc", 1}, {"Fucose", 1}}));
  TEST_EQUAL(aliases.size(), 5)
  TEST_REAL_SIMILAR(findIon(aliases, "glycan:diagnostic:Neu5Gc").getMZ(), 308.09761)
  TEST_REAL_SIMILAR(findIon(aliases, "glycan:diagnostic:Fuc").getMZ(), 147.06519)
  for (const auto& ion : ions)
  {
    TEST_EQUAL(ion.charge, 1)
  }
  TEST_TRUE(generator.getFragments(Composition {}).empty())
}
END_SECTION

START_SECTION((deoxyhexose diagnostic aliases preserve canonical and stereochemical identity))
{
  Generator::Options options;
  options.add_b_ions = options.add_y_ions = false;
  Generator generator(options);
  for (const auto& symbol : {"dHex", "d-Hex"})
  {
    const auto ions = generator.getFragments(comp({{symbol, 1}}));
    TEST_EQUAL(ions.size(), 1)
    const auto& ion = findIon(ions, "glycan:diagnostic:d-Hex");
    TEST_REAL_SIMILAR(ion.getMZ(), 147.06519)
    TEST_EQUAL(ion.composition.components.size(), 1)
    TEST_EQUAL(std::get<std::string>(ion.composition.components[0].first), "d-Hex")
    TEST_EQUAL(*MzPAF::parse(ion.getAnnotation()).named_compound, ion.name)
  }
  // Aliases describe the same component; Fuc additionally specifies stereochemistry.
  TEST_EQUAL(generator.getFragments(comp({{"dHex", 1}, {"d-Hex", 2}})).size(), 1)
  const auto mixed = generator.getFragments(comp({{"Fuc", 1}, {"dHex", 1}}));
  TEST_EQUAL(mixed.size(), 2)
  TEST_REAL_SIMILAR(findIon(mixed, "glycan:diagnostic:Fuc").getMZ(), findIon(mixed, "glycan:diagnostic:d-Hex").getMZ())
}
END_SECTION

START_SECTION((bounded composition B / Y / C / Z ions, charge ranges and formulas))
{
  Generator::Options options;
  options.add_diagnostic_ions = false;
  options.max_composition_size = 1;
  options.add_c_ions = options.add_z_ions = true;
  Generator generator(options);
  const auto ions = generator.getFragments(comp({{"Hex", 2}, {"HexNAc", 1}}));
  // Two retained compositions * four series * two charges, plus Y0 (Z0 has zero neutral mass).
  TEST_EQUAL(ions.size(), 18)
  const auto& b = findIon(ions, "glycan:B:composition:Hex1");
  const auto& c = findIon(ions, "glycan:C:composition:Hex1");
  const auto& y = findIon(ions, "glycan:Y:composition:Hex1");
  const auto& z = findIon(ions, "glycan:Z:composition:Hex1");
  TEST_REAL_SIMILAR(b.getMZ(), 163.06010)
  TEST_REAL_SIMILAR(c.getMZ() - b.getMZ(), 18.010565)
  TEST_REAL_SIMILAR(y.getMZ() - z.getMZ(), 18.010565)
  TEST_REAL_SIMILAR(findIon(ions, b.name, 2).getMZ(), (b.getMZ() + Constants::PROTON_MASS_U) / 2)
  TEST_FALSE(b.root_cleavage.has_value())
  TEST_FALSE(b.attachment_position.has_value())
  Composition formula;
  formula.components.emplace_back(ProForma::FormulaTag {"C6H10O5", std::nullopt}, 2);
  const auto custom = generator.getFragments(formula);
  TEST_REAL_SIMILAR(findIon(custom, "glycan:B:composition:Formula(C6H10O5)1").getMZ(), b.getMZ())
  const auto duplicate = generator.getFragments(comp({{"Hex", 1}, {"Hex", 1}, {"HexNAc", 1}, {"Hex", 0}}));
  TEST_EQUAL(duplicate.size(), ions.size())
  for (Size i = 0; i < ions.size(); ++i)
  {
    TEST_EQUAL(duplicate[i].name, ions[i].name)
    TEST_REAL_SIMILAR(duplicate[i].getMZ(), ions[i].getMZ())
  }
}
END_SECTION

START_SECTION((single neutral losses are feasible, independent and residue - specific))
{
  Generator::Options options;
  options.add_diagnostic_ions = options.add_y_ions = false;
  options.max_composition_size = 1;
  options.max_charge = 1;
  options.neutral_losses = {EmpiricalFormula("H2O"), EmpiricalFormula("C100"), EmpiricalFormula("H2O")};
  options.specific_neutral_losses["HexNAc"] = {EmpiricalFormula("C2H4O2")};
  const auto ions = Generator(options).getFragments(comp({{"Hex", 1}, {"HexNAc", 1}}));
  TEST_EQUAL(ions.size(), 5)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:B:composition:HexNAc1-H2O1").getMZ(), 186.07608)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:B:composition:HexNAc1-C2H4O2").getMZ(), 144.06552)
  for (const auto& ion : ions)
  {
    TEST_TRUE(ion.name.find("C100") == std::string::npos)
  }
}
END_SECTION

START_SECTION((branched tree fragments preserve cleavage identity and connectivity))
{
  GlycanStructure tree;
  tree.addMonosaccharide(std::string("HexNAc")); // 0
  tree.addMonosaccharide(std::string("Hex"), 0); // 1
  tree.addMonosaccharide(std::string("Hex"), 0); // 2
  tree.addMonosaccharide(std::string("Fuc"), 1); // 3
  Generator::Options options;
  options.add_diagnostic_ions = false;
  options.max_charge = 1;
  options.add_c_ions = options.add_z_ions = true;
  const auto ions = Generator(options).getFragments(tree);
  const auto& internal = findIon(ions, "glycan:B:tree:root=1;cuts=3;retained=Hex1");
  TEST_EQUAL(*internal.root_cleavage, 1)
  TEST_EQUAL(internal.branch_cleavages.size(), 1)
  TEST_EQUAL(internal.branch_cleavages[0], 3)
  TEST_REAL_SIMILAR(internal.getMZ(), 163.06010)
  const auto& y = findIon(ions, "glycan:Y:tree:root=0;cuts=1,2;retained=HexNAc1");
  TEST_REAL_SIMILAR(y.neutral_mass, EmpiricalFormula("C8H15NO6").getMonoWeight())
  const auto& z = findIon(ions, "glycan:Z:tree:root=0;cuts=1,2;retained=HexNAc1");
  TEST_REAL_SIMILAR(y.neutral_mass - z.neutral_mass, 2 * EmpiricalFormula("H2O").getMonoWeight())
  // Node 2 is a distinct, isobaric terminal fragment; do not collapse its interpretation.
  TEST_REAL_SIMILAR(findIon(ions, "glycan:B:tree:root=2;cuts=;retained=Hex1").getMZ(), internal.getMZ())
  options.add_internal_fragments = false;
  for (const auto& ion : Generator(options).getFragments(tree))
  {
    if (ion.ion_type == Ion::B || ion.ion_type == Ion::C) { TEST_TRUE(ion.branch_cleavages.empty()) }
  }
  options.allow_structural = false;
  const auto fallback = Generator(options).getFragments(tree);
  const auto composition = Generator(options).getFragments(tree.getComposition());
  TEST_EQUAL(fallback.size(), composition.size())
  for (Size i = 0; i < fallback.size(); ++i)
  {
    TEST_EQUAL(fallback[i].name, composition[i].name)
  }
  options.allow_structural = true;
  options.max_cleavages = 1;
  for (const auto& ion : Generator(options).getFragments(tree))
  {
    TEST_TRUE(ion.branch_cleavages.size() + (ion.root_cleavage.has_value() ? 1 : 0) <= 1)
  }
}
END_SECTION

START_SECTION((GlyPy stored fragment masses match the supported glycosidic series))
{
  // Adapted Apache-2.0 fixtures, pinned sources and notices are in the data directory.
  TOLERANCE_ABSOLUTE(0.0001)
  TOLERANCE_RELATIVE(1.0)
  Generator::Options options;
  options.add_diagnostic_ions = false;
  options.add_c_ions = options.add_z_ions = true;
  options.max_cleavages = 3;
  options.max_charge = 1;
  const auto cases = loadGlypyCases();
  std::map<std::string, std::vector<double>> observed, expected;
  const std::map<Ion, char> series {{Ion::B, 'B'}, {Ion::C, 'C'}, {Ion::Y, 'Y'}, {Ion::Z, 'Z'}};
  for (const auto& fragment : Generator(options).getFragments(cases.at("common_glycan").tree))
  {
    // GlyPy has no virtual bond below the reducing-end residue.
    if (fragment.root_cleavage == 0 || fragment.branch_cleavages == std::vector<Size> {0}) { continue; }
    const auto kind = fragment.root_cleavage ? std::string(1, series.at(fragment.ion_type)) + std::string(fragment.branch_cleavages.size(), 'Y')
                                             : std::string(fragment.branch_cleavages.size(), series.at(fragment.ion_type));
    observed[kind].push_back(fragment.neutral_mass);
  }
  std::ifstream input(OPENMS_GET_TEST_DATA_PATH("TheoreticalGlycanSpectrumGenerator_glypy/stored_masses.tsv"));
  TEST_TRUE(input.good())
  Size count = 0;
  for (std::string line; std::getline(input, line);)
  {
    if (line.empty() || line[0] == '#') { continue; }
    std::istringstream row(line);
    std::string kind;
    double mass;
    row >> kind >> mass;
    if (! row) { throw std::runtime_error("Malformed GlyPy stored mass: " + line); }
    expected[kind].push_back(mass);
    ++count;
  }
  TEST_EQUAL(count, 96)
  TEST_EQUAL(observed.size(), expected.size())
  for (auto& [kind, masses] : expected)
  {
    STATUS("GlyPy stored kind: " << kind)
    auto& actual = observed[kind];
    TEST_EQUAL(actual.size(), masses.size())
    std::sort(actual.begin(), actual.end());
    std::sort(masses.begin(), masses.end());
    for (Size i = 0; i < std::min(actual.size(), masses.size()); ++i)
    {
      TEST_REAL_SIMILAR(actual[i], masses[i])
    }
  }
  TOLERANCE_ABSOLUTE(0.00001)
  TOLERANCE_RELATIVE(1.00000001)
}
END_SECTION

START_SECTION((GlyPy tree fixtures preserve cleavage identities and masses at charges one to three))
{
  TOLERANCE_ABSOLUTE(0.0001)
  TOLERANCE_RELATIVE(1.0)
  const auto cases = loadGlypyCases();
  TEST_EQUAL(cases.size(), 3)
  Size comparisons = 0;
  for (const auto& [fixture, reference] : cases)
  {
    for (Size max_cleavages = 1; max_cleavages <= 3; ++max_cleavages)
    {
      STATUS("GlyPy fixture: " << fixture << "; maximum cleavages: " << max_cleavages)
      Generator::Options options;
      options.add_diagnostic_ions = false;
      options.add_c_ions = options.add_z_ions = true;
      options.max_charge = 3;
      options.max_cleavages = max_cleavages;
      using ChargedKey = std::pair<GlypyKey, Int>;
      std::map<ChargedKey, Generator::Fragment> observed;
      for (const auto& fragment : Generator(options).getFragments(reference.tree))
      {
        if (fragment.root_cleavage == 0 || fragment.branch_cleavages == std::vector<Size> {0}) { continue; }
        const GlypyKey key {fragment.ion_type, fragment.root_cleavage, fragment.branch_cleavages};
        const auto inserted = observed.emplace(ChargedKey {key, fragment.charge}, fragment);
        TEST_TRUE(inserted.second)
      }
      Size expected_count = 0;
      for (const auto& [key, mass] : reference.masses)
      {
        const auto& [series, root, branches] = key;
        if (branches.size() + root.has_value() > max_cleavages) { continue; }
        for (Int charge = 1; charge <= 3; ++charge)
        {
          ++expected_count;
          const auto found = observed.find(ChargedKey {key, charge});
          TEST_TRUE(found != observed.end())
          if (found == observed.end()) { continue; }
          TEST_REAL_SIMILAR(found->second.neutral_mass, mass)
          // Independent proton mass used by Pyteomics 5.0.1 (see fixture README).
          TEST_REAL_SIMILAR(found->second.getMZ(), mass / charge + 1.00727646677)
          ++comparisons;
        }
      }
      TEST_EQUAL(observed.size(), expected_count)
    }
  }
  TEST_EQUAL(comparisons, 2892)
  TOLERANCE_ABSOLUTE(0.00001)
  TOLERANCE_RELATIVE(1.00000001)
}
END_SECTION

START_SECTION((glycopeptide attachment, HCD / ETD / EThcD and configurable retention))
{
  const auto peptide = AASequence::fromString("ANST");
  const auto glycan = comp({{"HexNAc", 2}, {"Hex", 3}, {"Fuc", 1}});
  Generator generator;
  const auto hcd = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::HCD);
  const double stub_mass = EmpiricalFormula("C8H13NO5").getMonoWeight();
  const double full_mass = 2 * stub_mass + 3 * EmpiricalFormula("C6H10O5").getMonoWeight() + EmpiricalFormula("C6H10O4").getMonoWeight();
  const auto& stub = findIon(hcd, "peptide:b2;glycan=HexNAc1;site=N2");
  TEST_REAL_SIMILAR(stub.neutral_mass, peptide.getPrefix(2).getMonoWeight(Residue::BIon) + stub_mass)
  TEST_EQUAL(*stub.attachment_position, 1)
  TEST_EQUAL(stub.attachment_residue, "N")
  TEST_REAL_SIMILAR(findIon(hcd, "peptide:b2;glycan=0;site=N2").neutral_mass, peptide.getPrefix(2).getMonoWeight(Residue::BIon))
  TEST_FALSE(findIon(hcd, "peptide:b1").attachment_position.has_value())
  TEST_FALSE(findIon(hcd, "peptide:y2").attachment_position.has_value())
  TEST_REAL_SIMILAR(findIon(hcd, "glycan:Y:composition:0;site=N2").neutral_mass, peptide.getMonoWeight())
  const auto etd = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::ETD);
  TEST_EQUAL(etd.size(), 12)
  const auto& intact = findIon(etd, "peptide:c2;glycan=Fuc1Hex3HexNAc2;site=N2");
  TEST_REAL_SIMILAR(intact.neutral_mass, peptide.getPrefix(2).getMonoWeight(Residue::CIon) + full_mass)
  // The radical NST suffix has neutral formula C11H18N3O7.
  TEST_REAL_SIMILAR(findIon(etd, "peptide:z3;glycan=Fuc1Hex3HexNAc2;site=N2").neutral_mass, 304.11447493 + full_mass)
  for (const auto& ion : etd)
  {
    TEST_TRUE(ion.ion_type == Ion::PEPTIDE)
  }
  const auto ethcd = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::ETHCD);
  const Size glycan_ions = std::count_if(hcd.begin(), hcd.end(), [](const auto& ion) { return ion.ion_type != Ion::PEPTIDE; });
  TEST_EQUAL(ethcd.size(), glycan_ions + 2 * etd.size())
  TEST_REAL_SIMILAR(findIon(ethcd, intact.name).getMZ(), intact.getMZ())
  TEST_REAL_SIMILAR(findIon(ethcd, "peptide:b2;glycan=Fuc1Hex3HexNAc2;site=N2").neutral_mass,
                    peptide.getPrefix(2).getMonoWeight(Residue::BIon) + full_mass)
  auto options = generator.getOptions();
  Generator::PeptideRetention rule;
  rule.stripped = false;
  rule.stubs = {comp({{"HexNAc", 1}, {"Fuc", 1}})};
  options.peptide_retention['b'] = rule;
  generator.setOptions(options);
  const auto custom = generator.getGlycopeptideFragments(peptide, glycan, 1, Generator::FragmentationMethod::HCD);
  TEST_REAL_SIMILAR(findIon(custom, "peptide:b2;glycan=Fuc1HexNAc1;site=N2").neutral_mass,
                    peptide.getPrefix(2).getMonoWeight(Residue::BIon) + stub_mass + EmpiricalFormula("C6H10O4").getMonoWeight())
  // No available fucose: explicit retention rules must fail instead of inventing mass.
  TEST_EXCEPTION(Exception::InvalidParameter,
                 generator.getGlycopeptideFragments(peptide, comp({{"HexNAc", 1}}), 1, Generator::FragmentationMethod::HCD))
  TEST_EXCEPTION(Exception::InvalidParameter, generator.getGlycopeptideFragments(peptide, glycan, 4, Generator::FragmentationMethod::HCD))
}
END_SECTION

START_SECTION((ETD and EThcD radical backbone masses match independent elemental expectations))
{
  Generator generator;
  const auto peptide = AASequence::fromString("ANST");
  const auto glycan = comp({{"HexNAc", 1}});
  // Fixed monoisotopic masses from elemental counts, without AASequence's ion conversions:
  // z1(T): C4H7O3; z3(NST)+HexNAc: C19H31N4O12; c2(AN)+HexNAc: C15H27N5O8.
  // Radical z = internal residue sum + O - N; charge adds one proton per unit.
  for (const auto method : {Generator::FragmentationMethod::ETD, Generator::FragmentationMethod::ETHCD})
  {
    const auto ions = generator.getGlycopeptideFragments(peptide, glycan, 1, method);
    const auto& z1 = findIon(ions, "peptide:z1");
    TEST_FALSE(z1.attachment_position.has_value())
    TEST_REAL_SIMILAR(z1.neutral_mass, 103.03951908)
    TEST_REAL_SIMILAR(z1.getMZ(), 104.04679555)
    TEST_REAL_SIMILAR(findIon(ions, "peptide:z1", 2).getMZ(), 52.52703601)
    const auto& z3 = findIon(ions, "peptide:z3;glycan=HexNAc1;site=N2");
    TEST_EQUAL(*z3.attachment_position, 1)
    TEST_REAL_SIMILAR(z3.neutral_mass, 507.19384745)
    TEST_REAL_SIMILAR(z3.getMZ(), 508.20112392)
    TEST_REAL_SIMILAR(findIon(ions, z3.name, 2).getMZ(), 254.60420019)
    TEST_REAL_SIMILAR(findIon(ions, "peptide:c2;glycan=HexNAc1;site=N2").getMZ(), 406.19323932)
    auto modified = peptide;
    modified.setCTerminalModification("Amidated"); // OH -> NH2: -0.98401558 Da.
    const auto amidated = generator.getGlycopeptideFragments(modified, glycan, 1, method);
    TEST_REAL_SIMILAR(findIon(amidated, "peptide:z1").getMZ(), 103.06277997)
    TEST_REAL_SIMILAR(findIon(amidated, z3.name).getMZ(), 507.21710834)
    TEST_REAL_SIMILAR(findIon(amidated, z3.name, 2).getMZ(), 254.11219240)
    TEST_REAL_SIMILAR(findIon(amidated, "peptide:c2;glycan=HexNAc1;site=N2").getMZ(), 406.19323932)
  }
}
END_SECTION

START_SECTION((fragment annotations reject names that cannot roundtrip through mzPAF))
{
  Generator::Fragment fragment;
  for (const auto& name :
       std::vector<std::string> {"", "glycan name", "glycan\tname", "glycan\nname", "glycan[name", "glycan]name", std::string("glycan\0name", 11)})
  {
    fragment.name = name;
    TEST_EXCEPTION(Exception::InvalidParameter, fragment.getAnnotation())
  }
  fragment.name = "glycan:diagnostic:d-Hex";
  fragment.charge = 2;
  fragment.neutral_mass = 146.0579088;
  const auto parsed = MzPAF::parse(fragment.getAnnotation());
  TEST_EQUAL(*parsed.named_compound, fragment.name)
  TEST_EQUAL(*parsed.charge, 2)
  // The public C++ fragment type can also reach serialization through spectrum export.
  fragment.name = "glycan name";
  TEST_EXCEPTION(Exception::InvalidParameter, Generator::toSpectrum({fragment}))
}
END_SECTION

START_SECTION((spectrum metadata and mzPAF roundtrip))
{
  auto ions = Generator().getFragments(comp({{"HexNAc", 2}, {"Hex", 1}}));
  std::reverse(ions.begin(), ions.end());
  const auto spectrum = Generator::toSpectrum(ions);
  TEST_EQUAL(spectrum.getMSLevel(), 2)
  TEST_EQUAL(spectrum.size(), ions.size())
  TEST_TRUE(spectrum.isSorted())
  TEST_EQUAL(spectrum.getStringDataArrays()[0].getName(), Constants::UserParam::IonNames)
  TEST_EQUAL(spectrum.getIntegerDataArrays()[0].getName(), "Charges")
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    const auto annotation = MzPAF::parse(spectrum.getStringDataArrays()[0][i]);
    TEST_TRUE(annotation.ion_series == MzPAFIonSeries::NAMED)
    TEST_EQUAL(annotation.charge.value_or(1), spectrum.getIntegerDataArrays()[0][i])
    const auto& ion = findIon(ions, *annotation.named_compound, annotation.charge.value_or(1));
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), ion.getMZ())
  }
}
END_SECTION

START_SECTION((validation and bounded work))
{
  Generator generator;
  TEST_EXCEPTION(Exception::InvalidParameter, generator.getFragments(comp({{"Hex", -1}})))
  TEST_EXCEPTION(Exception::ElementNotFound, generator.getFragments(comp({{"Unknown", 1}})))
  Generator::Options options;
  options.max_charge = 0;
  TEST_EXCEPTION(Exception::InvalidParameter, generator.setOptions(options))
  TEST_EQUAL(generator.getOptions().max_charge, 2)
  options = Generator::Options {};
  options.max_fragments = 1;
  TEST_EXCEPTION(Exception::InvalidParameter, Generator(options).getFragments(comp({{"HexNAc", 1}})))
  options = Generator::Options {};
  options.max_states = 1;
  TEST_EXCEPTION(Exception::InvalidParameter, Generator(options).getFragments(comp({{"Hex", 1000000}})))
  options = Generator::Options {};
  options.add_b_ions = options.add_y_ions = false;
  options.max_oxonium_charge = 2;
  const auto ions = Generator(options).getFragments(comp({{"HexNAc", 1000000}}));
  TEST_EQUAL(ions.size(), 12)
  TEST_REAL_SIMILAR(findIon(ions, "glycan:diagnostic:HexNAc", 2).getMZ(), (204.08665 + Constants::PROTON_MASS_U) / 2)
}
END_SECTION

END_TEST
