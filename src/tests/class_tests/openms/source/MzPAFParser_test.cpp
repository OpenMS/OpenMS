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

#include <OpenMS/CHEMISTRY/MzPAFParser.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideHit.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MzPAFParser, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(Simple ion parsing - y ion)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.has_value(), true)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.charge.has_value(), false)
  TEST_EQUAL(ann.isValid(), true)
}
END_SECTION

START_SECTION(Simple ion parsing - b ion)
{
  MzPAFAnnotation ann = MzPAFParser::parse("b2");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.isValid(), true)
}
END_SECTION

START_SECTION(Simple ion parsing - all standard ions)
{
  // a-ion
  MzPAFAnnotation a = MzPAFParser::parse("a3");
  TEST_EQUAL(a.ion_series, MzPAFIonSeries::A)
  TEST_EQUAL(a.ordinal.value(), 3)

  // c-ion
  MzPAFAnnotation c = MzPAFParser::parse("c5");
  TEST_EQUAL(c.ion_series, MzPAFIonSeries::C)
  TEST_EQUAL(c.ordinal.value(), 5)

  // x-ion
  MzPAFAnnotation x = MzPAFParser::parse("x1");
  TEST_EQUAL(x.ion_series, MzPAFIonSeries::X)
  TEST_EQUAL(x.ordinal.value(), 1)

  // z-ion
  MzPAFAnnotation z = MzPAFParser::parse("z7");
  TEST_EQUAL(z.ion_series, MzPAFIonSeries::Z)
  TEST_EQUAL(z.ordinal.value(), 7)
}
END_SECTION

START_SECTION(Ion with charge)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4^2");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.charge.has_value(), true)
  TEST_EQUAL(ann.charge.value(), 2)
}
END_SECTION

START_SECTION(Ion with neutral loss)
{
  MzPAFAnnotation ann = MzPAFParser::parse("b2-H2O");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.neutral_losses.size(), 1)
  TEST_EQUAL(ann.neutral_losses[0].original_text, "H2O")
}
END_SECTION

START_SECTION(Ion with multiple neutral losses)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y5-H2O-NH3");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 5)
  TEST_EQUAL(ann.neutral_losses.size(), 2)
  TEST_EQUAL(ann.neutral_losses[0].original_text, "H2O")
  TEST_EQUAL(ann.neutral_losses[1].original_text, "NH3")
}
END_SECTION

START_SECTION(Ion with isotope offset)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y2+2i");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.isotope_offset.has_value(), true)
  TEST_EQUAL(ann.isotope_offset.value(), 2)
}
END_SECTION

START_SECTION(Ion with mass delta in Da)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4/0.001");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.mass_delta.has_value(), true)
  TEST_REAL_SIMILAR(ann.mass_delta.value().value, 0.001)
  TEST_EQUAL(ann.mass_delta.value().unit, MzPAFDeltaUnit::DALTON)
}
END_SECTION

START_SECTION(Ion with mass delta in ppm)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4/-1.4ppm");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.mass_delta.has_value(), true)
  TEST_REAL_SIMILAR(ann.mass_delta.value().value, -1.4)
  TEST_EQUAL(ann.mass_delta.value().unit, MzPAFDeltaUnit::PPM)
}
END_SECTION

START_SECTION(Ion with confidence)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4*0.75");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 4)
  TEST_EQUAL(ann.confidence.has_value(), true)
  TEST_REAL_SIMILAR(ann.confidence.value(), 0.75)
}
END_SECTION

START_SECTION(Complex annotation with all modifiers)
{
  MzPAFAnnotation ann = MzPAFParser::parse("b2-H2O^2/3.2ppm*0.75");
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
  MzPAFAnnotation ann = MzPAFParser::parse("IY");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::IMMONIUM)
  TEST_EQUAL(ann.immonium_residue.has_value(), true)
  TEST_EQUAL(ann.immonium_residue.value(), 'Y')
}
END_SECTION

START_SECTION(Internal fragment)
{
  MzPAFAnnotation ann = MzPAFParser::parse("m3:6");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::INTERNAL)
  TEST_EQUAL(ann.internal_range.has_value(), true)
  TEST_EQUAL(ann.internal_range.value().first, 3)
  TEST_EQUAL(ann.internal_range.value().second, 6)
}
END_SECTION

START_SECTION(Reporter ion)
{
  MzPAFAnnotation ann = MzPAFParser::parse("r[TMT127N]");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::REPORTER)
  TEST_EQUAL(ann.reporter_name.has_value(), true)
  TEST_EQUAL(ann.reporter_name.value(), "TMT127N")
}
END_SECTION

START_SECTION(Formula ion)
{
  MzPAFAnnotation ann = MzPAFParser::parse("f{C16H22O}");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::FORMULA)
  TEST_EQUAL(ann.formula.has_value(), true)
  // EmpiricalFormula.toString() includes count even if 1 (C16H22O1)
  TEST_EQUAL(ann.formula.value().toString(), "C16H22O1")
}
END_SECTION

START_SECTION(Precursor ion)
{
  MzPAFAnnotation ann = MzPAFParser::parse("p");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::PRECURSOR)
}
END_SECTION

START_SECTION(Multi-analyte annotation)
{
  MzPAFAnnotation ann = MzPAFParser::parse("1@y12");
  TEST_EQUAL(ann.analyte_index.has_value(), true)
  TEST_EQUAL(ann.analyte_index.value(), 1)
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(ann.ordinal.value(), 12)
}
END_SECTION

START_SECTION(Embedded sequence)
{
  MzPAFAnnotation ann = MzPAFParser::parse("b2{LC}");
  TEST_EQUAL(ann.ion_series, MzPAFIonSeries::B)
  TEST_EQUAL(ann.ordinal.value(), 2)
  TEST_EQUAL(ann.embedded_sequence.has_value(), true)
  TEST_EQUAL(ann.embedded_sequence.value(), "LC")
}
END_SECTION

START_SECTION(Parse multiple annotations)
{
  MzPAFPeakAnnotations anns = MzPAFParser::parseMultiple("b2,y4^2");
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
  auto valid = MzPAFParser::tryParse("y4");
  TEST_EQUAL(valid.has_value(), true)

  auto invalid = MzPAFParser::tryParse("invalid[");
  TEST_EQUAL(invalid.has_value(), false)
}
END_SECTION

START_SECTION(toString LOSSLESS mode)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4^2-H2O/0.001*0.75");
  String s = MzPAFParser::toString(ann, MzPAFWriteMode::LOSSLESS);
  TEST_EQUAL(s, "y4-H2O^2/0.001*0.75")
}
END_SECTION

START_SECTION(toString CANONICAL mode)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4^2");
  String s = MzPAFParser::toString(ann, MzPAFWriteMode::CANONICAL);
  TEST_EQUAL(s, "y4^2")
}
END_SECTION

START_SECTION(toString multiple annotations)
{
  MzPAFPeakAnnotations anns = MzPAFParser::parseMultiple("b2,y4^2");
  String s = MzPAFParser::toString(anns, MzPAFWriteMode::CANONICAL);
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
    auto ann = MzPAFParser::parse(input);
    String output = MzPAFParser::toString(ann, MzPAFWriteMode::CANONICAL);
    auto reparsed = MzPAFParser::parse(output);
    TEST_EQUAL(ann, reparsed)
  }
}
END_SECTION

START_SECTION(isMzPAFFormat detection)
{
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("y4"), true)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("b2-H2O"), true)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("y4^2"), true)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("IY"), true)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("m3:6"), true)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("r[TMT127N]"), true)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat(""), false)
  TEST_EQUAL(MzPAFParser::isMzPAFFormat("random text"), false)
}
END_SECTION

START_SECTION(toPeakAnnotation integration)
{
  MzPAFAnnotation ann = MzPAFParser::parse("y4^2");
  PeptideHit::PeakAnnotation pa = MzPAFParser::toPeakAnnotation(ann, 500.123, 1000.0);

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

  MzPAFPeakAnnotations anns = MzPAFParser::fromPeakAnnotation(pa);
  TEST_EQUAL(anns.size(), 1)
  TEST_EQUAL(anns.annotations[0].ion_series, MzPAFIonSeries::Y)
  TEST_EQUAL(anns.annotations[0].ordinal.value(), 4)
  TEST_EQUAL(anns.annotations[0].charge.value(), 2)
  TEST_EQUAL(anns.annotations[0].neutral_losses.size(), 1)
}
END_SECTION

START_SECTION(ionSeriesToChar and charToIonSeries)
{
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::A), 'a')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::B), 'b')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::C), 'c')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::X), 'x')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::Y), 'y')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::Z), 'z')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::PRECURSOR), 'p')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::IMMONIUM), 'I')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::INTERNAL), 'm')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::REPORTER), 'r')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::FORMULA), 'f')
  TEST_EQUAL(MzPAFParser::ionSeriesToChar(MzPAFIonSeries::NAMED), '_')

  MzPAFIonSeries series;
  TEST_EQUAL(MzPAFParser::charToIonSeries('y', series), true)
  TEST_EQUAL(series, MzPAFIonSeries::Y)

  TEST_EQUAL(MzPAFParser::charToIonSeries('Q', series), false)
}
END_SECTION

START_SECTION(Error handling - empty input)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAFParser::parse(""))
}
END_SECTION

START_SECTION(Error handling - invalid ion series)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAFParser::parse("Q5"))
}
END_SECTION

START_SECTION(Error handling - unclosed bracket)
{
  TEST_EXCEPTION(MzPAFParseError, MzPAFParser::parse("r[TMT127N"))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
