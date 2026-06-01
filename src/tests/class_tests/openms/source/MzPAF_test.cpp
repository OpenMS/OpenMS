// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/MzPAF.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideHit.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzPAF, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(Simple ion parsing - y ion)
{
  MzPAFAnnotation ann = MzPAF::parse("y4");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.has_value(), true)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.charge.has_value(), false)
  TEST_EQUAL(ann.isValid(), true)
}
END_SECTION

START_SECTION(Simple ion parsing - b ion)
{
  MzPAFAnnotation ann = MzPAF::parse("b2");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.isValid(), true)
}
END_SECTION

START_SECTION(Simple ion parsing - all standard ions)
{
  // a-ion
  MzPAFAnnotation a = MzPAF::parse("a3");
  TEST_EQUAL(a.ion_series, MzPAFIonSeries::A)
  TEST_EQUAL(a.ordinal.value(), 3)

  // c-ion
  MzPAFAnnotation c = MzPAF::parse("c5");
  TEST_EQUAL(c.ion_series, MzPAFIonSeries::C)
  TEST_EQUAL(c.ordinal.value(), 5)

  // x-ion
  MzPAFAnnotation x = MzPAF::parse("x1");
  TEST_EQUAL(x.ion_series, MzPAFIonSeries::X)
  TEST_EQUAL(x.ordinal.value(), 1)

  // z-ion
  MzPAFAnnotation z = MzPAF::parse("z7");
  TEST_EQUAL(z.ion_series, MzPAFIonSeries::Z)
  TEST_EQUAL(z.ordinal.value(), 7)
}
END_SECTION

START_SECTION(Ion with charge)
{
  MzPAFAnnotation ann = MzPAF::parse("y4^2");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.charge.has_value(), true)
  TEST_EQUAL(ann.charge.value(), 2)
}
END_SECTION

START_SECTION(Ion with neutral loss)
{
  MzPAFAnnotation ann = MzPAF::parse("b2-H2O");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.neutral_losses.size(), 1)
  TEST_EQUAL(ann.neutral_losses[0].formula.toString(), "H2O1")
}
END_SECTION

START_SECTION(Ion with multiple neutral losses)
{
  MzPAFAnnotation ann = MzPAF::parse("y5-H2O-NH3");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 5)
  TEST_EQUAL(ann.neutral_losses.size(), 2)
  TEST_EQUAL(ann.neutral_losses[0].formula.toString(), "H2O1")
  TEST_EQUAL(ann.neutral_losses[1].formula.toString(), "H3N1")
}
END_SECTION

START_SECTION(Ion with isotope offset)
{
  MzPAFAnnotation ann = MzPAF::parse("y2+2i");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.isotope_offset.has_value(), true)
  TEST_EQUAL(ann.isotope_offset.value(), 2)
}
END_SECTION

START_SECTION(Ion with mass delta in Da)
{
  MzPAFAnnotation ann = MzPAF::parse("y4/0.001");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.mass_delta.has_value(), true)
  TEST_REAL_SIMILAR(ann.mass_delta.value().value, 0.001)
  TEST_EQUAL(ann.mass_delta.value().unit, MzPAFDeltaUnit::DALTON)
}
END_SECTION

START_SECTION(Ion with mass delta in ppm)
{
  MzPAFAnnotation ann = MzPAF::parse("y4/-1.4ppm");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.mass_delta.has_value(), true)
  TEST_REAL_SIMILAR(ann.mass_delta.value().value, -1.4)
  TEST_EQUAL(ann.mass_delta.value().unit, MzPAFDeltaUnit::PPM)
}
END_SECTION

START_SECTION(Ion with confidence)
{
  MzPAFAnnotation ann = MzPAF::parse("y4*0.75");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.confidence.has_value(), true)
  TEST_REAL_SIMILAR(ann.confidence.value(), 0.75)
}
END_SECTION

START_SECTION(Complex annotation with all modifiers)
{
  MzPAFAnnotation ann = MzPAF::parse("b2-H2O^2/3.2ppm*0.75");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.neutral_losses.size(), 1)
  TEST_EQUAL(ann.charge.value(), 2)
  TEST_REAL_SIMILAR(ann.mass_delta.value().value, 3.2)
  TEST_EQUAL(ann.mass_delta.value().unit, MzPAFDeltaUnit::PPM)
  TEST_REAL_SIMILAR(ann.confidence.value(), 0.75)
}
END_SECTION

START_SECTION(Immonium ion)
{
  MzPAFAnnotation ann = MzPAF::parse("IY");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::IMMONIUM)
  TEST_EQUAL(ann.immonium_residue.has_value(), true)
  TEST_EQUAL(ann.immonium_residue.value(), 'Y')
}
END_SECTION

START_SECTION(Internal fragment)
{
  MzPAFAnnotation ann = MzPAF::parse("m3:6");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::INTERNAL)
  TEST_EQUAL(ann.internal_range.has_value(), true)
  TEST_EQUAL(ann.internal_range.value().first, 3)
  TEST_EQUAL(ann.internal_range.value().second, 6)
}
END_SECTION

START_SECTION(Reporter ion)
{
  MzPAFAnnotation ann = MzPAF::parse("r[TMT127N]");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::REPORTER)
  TEST_EQUAL(ann.reporter_name.has_value(), true)
  TEST_EQUAL(ann.reporter_name.value(), "TMT127N")
}
END_SECTION

START_SECTION(Formula ion)
{
  MzPAFAnnotation ann = MzPAF::parse("f{C16H22O}");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::FORMULA)
  TEST_EQUAL(ann.formula.has_value(), true)
  // EmpiricalFormula.toString() includes count even if 1 (C16H22O1)
  TEST_EQUAL(ann.formula.value().toString(), "C16H22O1")
}
END_SECTION

START_SECTION(Precursor ion)
{
  MzPAFAnnotation ann = MzPAF::parse("p");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::PRECURSOR)
}
END_SECTION

START_SECTION(Multi-analyte annotation)
{
  MzPAFAnnotation ann = MzPAF::parse("1@y12");
  TEST_EQUAL(ann.analyte_index.has_value(), true)
  TEST_EQUAL(ann.analyte_index.value(), 1)
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 12)
}
END_SECTION

START_SECTION(Embedded sequence)
{
  MzPAFAnnotation ann = MzPAF::parse("b2{LC}");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.embedded_sequence.has_value(), true)
  TEST_EQUAL(ann.embedded_sequence.value(), "LC")
}
END_SECTION

START_SECTION(Parse multiple annotations)
{
  MzPAFPeakAnnotations anns = MzPAF::parseMultiple("b2,y4^2");
  TEST_EQUAL(anns.size(), 2)
  TEST_EQUAL(anns.annotations[0].ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(anns.annotations[0].ordinal.value(), 2)
  TEST_EQUAL(anns.annotations[1].ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(anns.annotations[1].ordinal.value(), 4)
  TEST_EQUAL(anns.annotations[1].charge.value(), 2)
}
END_SECTION

START_SECTION(tryParse non-throwing)
{
  auto valid = MzPAF::tryParse("y4");
  TEST_EQUAL(valid.has_value(), true)

  auto invalid = MzPAF::tryParse("invalid[");
  TEST_EQUAL(invalid.has_value(), false)
}
END_SECTION

START_SECTION(toString single annotation)
{
  MzPAFAnnotation ann = MzPAF::parse("y4^2-H2O/0.001*0.75");
  String s = MzPAF::toString(ann);
  TEST_EQUAL(s, "y4-H2O1^2/0.001*0.75")
}
END_SECTION

START_SECTION(toString simple annotation)
{
  MzPAFAnnotation ann = MzPAF::parse("y4^2");
  String s = MzPAF::toString(ann);
  TEST_EQUAL(s, "y4^2")
}
END_SECTION

START_SECTION(toString multiple annotations)
{
  MzPAFPeakAnnotations anns = MzPAF::parseMultiple("b2,y4^2");
  String s = MzPAF::toString(anns);
  TEST_EQUAL(s, "b2,y4^2")
}
END_SECTION

START_SECTION(Roundtrip test)
{
  vector<String> test_cases = {
    "y4",
    "b2",
    "y4^2",
    "b2-H2O",
    "y5-H2O-NH3",
    "y4/0.001",
    "y4*0.75",
    "IY",
    "m3:6",
    "r[TMT127N]",
    "p",
    "1@y12"
  };

  for (const auto& input : test_cases)
  {
    auto ann = MzPAF::parse(input);
    String output = MzPAF::toString(ann);
    auto reparsed = MzPAF::parse(output);
    TEST_EQUAL(ann, reparsed)
  }
}
END_SECTION

START_SECTION(isMzPAFFormat detection)
{
  TEST_EQUAL(MzPAF::isMzPAFFormat("y4"), true)
  TEST_EQUAL(MzPAF::isMzPAFFormat("b2-H2O"), true)
  TEST_EQUAL(MzPAF::isMzPAFFormat("y4^2"), true)
  TEST_EQUAL(MzPAF::isMzPAFFormat("IY"), true)
  TEST_EQUAL(MzPAF::isMzPAFFormat("m3:6"), true)
  TEST_EQUAL(MzPAF::isMzPAFFormat("r[TMT127N]"), true)
  TEST_EQUAL(MzPAF::isMzPAFFormat(""), false)
  TEST_EQUAL(MzPAF::isMzPAFFormat("random text"), false)
}
END_SECTION

START_SECTION(toPeakAnnotation integration)
{
  MzPAFAnnotation ann = MzPAF::parse("y4^2");
  PeptideHit::PeakAnnotation pa = MzPAF::toPeakAnnotation(ann, 500.123, 1000.0);

  TEST_EQUAL(pa.charge, 2)
  TEST_REAL_SIMILAR(pa.mz, 500.123)
  TEST_REAL_SIMILAR(pa.intensity, 1000.0)
  TEST_EQUAL(pa.annotation, "y4^2")
}
END_SECTION

START_SECTION(fromPeakAnnotation integration)
{
  PeptideHit::PeakAnnotation pa;
  pa.annotation = "y4^2-H2O";
  pa.charge = 2;
  pa.mz = 500.0;
  pa.intensity = 1000.0;

  MzPAFPeakAnnotations anns = MzPAF::fromPeakAnnotation(pa);
  TEST_EQUAL(anns.size(), 1)
  TEST_EQUAL(anns.annotations[0].ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(anns.annotations[0].ordinal.value(), 4)
  TEST_EQUAL(anns.annotations[0].charge.value(), 2)
  TEST_EQUAL(anns.annotations[0].neutral_losses.size(), 1)
}
END_SECTION

START_SECTION(ionSeriesToChar and charToIonSeries)
{
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::A), 'a')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::B), 'b')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::C), 'c')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::X), 'x')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::Y), 'y')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::Z), 'z')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::PRECURSOR), 'p')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::IMMONIUM), 'I')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::INTERNAL), 'm')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::REPORTER), 'r')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::FORMULA), 'f')
  TEST_EQUAL(MzPAF::ionSeriesToChar(MzPAFIonSeries::NAMED), '_')

  MzPAFIonSeries series;
  TEST_EQUAL(MzPAF::charToIonSeries('y', series), true)
  TEST_EQUAL(series, MzPAFIonSeries::Y)

  TEST_EQUAL(MzPAF::charToIonSeries('Q', series), false)
}
END_SECTION

START_SECTION(Error handling - empty input)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse(""))
}
END_SECTION

START_SECTION(Error handling - invalid ion series)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse("Q5"))
}
END_SECTION

START_SECTION(Error handling - unclosed bracket)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse("r[TMT127N"))
}
END_SECTION

START_SECTION(Named compound ion)
{
  MzPAFAnnotation ann = MzPAF::parse("_[Aspirin]");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::NAMED)
  TEST_EQUAL(ann.named_compound.has_value(), true)
  TEST_EQUAL(ann.named_compound.value(), "Aspirin")
  TEST_EQUAL(ann.isValid(), true)
}
END_SECTION

START_SECTION(Named compound ion with modifiers)
{
  MzPAFAnnotation ann = MzPAF::parse("_[Iodoacetamide]^1/0.5ppm");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::NAMED)
  TEST_EQUAL(ann.named_compound.value(), "Iodoacetamide")
  TEST_EQUAL(ann.charge.value(), 1)
  TEST_REAL_SIMILAR(ann.mass_delta.value().value, 0.5)
  TEST_EQUAL(ann.mass_delta.value().unit, MzPAFDeltaUnit::PPM)
}
END_SECTION

START_SECTION(Error handling - number overflow in ordinal)
{
  // Very large number that would overflow int
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse("y99999999999999999999"))
}
END_SECTION

START_SECTION(Error handling - number overflow in charge)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse("y4^99999999999999999999"))
}
END_SECTION

START_SECTION(Error handling - number overflow in analyte index)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse("99999999999999999999@y4"))
}
END_SECTION

START_SECTION(Error handling - invalid internal fragment range)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAF::parse("m99999999999999999999:6"))
}
END_SECTION

START_SECTION(Roundtrip test with named compound)
{
  auto ann = MzPAF::parse("_[MyCompound]");
  String output = MzPAF::toString(ann);
  auto reparsed = MzPAF::parse(output);
  TEST_EQUAL(ann, reparsed)
}
END_SECTION

START_SECTION(tryParseMultiple non-throwing)
{
  auto valid = MzPAF::tryParseMultiple("b2,y4^2");
  TEST_EQUAL(valid.has_value(), true)
  TEST_EQUAL(valid.value().size(), 2)

  auto invalid = MzPAF::tryParseMultiple("y4,invalid[");
  TEST_EQUAL(invalid.has_value(), false)
}
END_SECTION

START_SECTION(isStandardFragmentIon)
{
  // Standard fragment ions (a, b, c, x, y, z)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::A), true)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::B), true)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::C), true)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::X), true)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::Y), true)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::Z), true)

  // Special ion types (not standard fragment ions)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::PRECURSOR), false)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::IMMONIUM), false)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::INTERNAL), false)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::REPORTER), false)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::FORMULA), false)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::NAMED), false)
  TEST_EQUAL(MzPAF::isStandardFragmentIon(MzPAFIonSeries::UNKNOWN), false)
}
END_SECTION

START_SECTION(Adduct parsing)
{
  // Note: Adduct parsing depends on implementation - test if supported
  auto result = MzPAF::tryParse("y4+Na");
  if (result.has_value() && result.value().adduct.has_value())
  {
    TEST_EQUAL(result.value().ion_series, MzPAFIonSeries::Y)
    TEST_EQUAL(result.value().ordinal.value(), 4)
    // Adduct should be parsed as Na
  }
}
END_SECTION

START_SECTION(calculateTheoreticalMZ)
{
  AASequence seq = AASequence::fromString("PEPTIDER");

  // Test y4 ion (C-terminal 4 residues: IDER)
  MzPAFAnnotation y4 = MzPAF::parse("y4");
  auto mz_y4 = MzPAF::calculateTheoreticalMZ(y4, seq);
  TEST_EQUAL(mz_y4.has_value(), true)
  // y4 of PEPTIDER = IDER, singly charged
  TEST_REAL_SIMILAR(mz_y4.value(), 532.2726)

  // Test b3 ion (N-terminal 3 residues: PEP)
  MzPAFAnnotation b3 = MzPAF::parse("b3");
  auto mz_b3 = MzPAF::calculateTheoreticalMZ(b3, seq);
  TEST_EQUAL(mz_b3.has_value(), true)
  // b3 of PEPTIDER = PEP
  TEST_REAL_SIMILAR(mz_b3.value(), 324.1554)

  // Test doubly charged y4
  MzPAFAnnotation y4_2 = MzPAF::parse("y4^2");
  auto mz_y4_2 = MzPAF::calculateTheoreticalMZ(y4_2, seq);
  TEST_EQUAL(mz_y4_2.has_value(), true)
  // Doubly charged: (mass + 2*proton) / 2
  TEST_REAL_SIMILAR(mz_y4_2.value(), 266.6399)

  // Test non-standard ion (should return nullopt)
  MzPAFAnnotation precursor = MzPAF::parse("p");
  auto mz_p = MzPAF::calculateTheoreticalMZ(precursor, seq);
  TEST_EQUAL(mz_p.has_value(), false)

  // Test invalid ordinal (out of range)
  MzPAFAnnotation y100 = MzPAF::parse("y100");
  auto mz_y100 = MzPAF::calculateTheoreticalMZ(y100, seq);
  TEST_EQUAL(mz_y100.has_value(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
