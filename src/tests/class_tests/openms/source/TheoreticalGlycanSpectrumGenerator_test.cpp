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
#include <algorithm>
#include <set>
#include <string_view>

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
  TEST_REAL_SIMILAR(findIon(etd, "peptide:z3;glycan=Fuc1Hex3HexNAc2;site=N2").neutral_mass,
                    peptide.getSuffix(3).getMonoWeight(Residue::ZIon) + full_mass)
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
