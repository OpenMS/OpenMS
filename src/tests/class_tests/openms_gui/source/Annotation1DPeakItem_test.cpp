// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>

#include <OpenMS/KERNEL/Peak1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

/// A label as it sits on the plot, converted back into the annotation stored in the PeptideHit.
/// This is the direction TOPPView takes whenever it writes edited annotations back to the ID data.
PeptideHit::PeakAnnotation toAnnotation_(const QString& label, int charge = 0)
{
  Annotation1DPeakItem<Peak1D> item(Peak1D(1234.5, 678.0), label, QColor(Qt::black), charge);
  return item.toPeakAnnotation();
}

START_TEST(Annotation1DPeakItem, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((PeptideHit::PeakAnnotation toPeakAnnotation() const))
{
  // the ion name keeps its charge notation and 'charge' only repeats it as a number
  // (see PeptideHit::PeakAnnotation). Stripping it here made the identification view append a second
  // suffix on the next redraw, which showed every fragment as 2+ (issue #8766).
  PeptideHit::PeakAnnotation fa = toAnnotation_("y5++");
  TEST_STRING_EQUAL(fa.annotation, "y5++")
  TEST_EQUAL(fa.charge, 2)
  TEST_REAL_SIMILAR(fa.mz, 1234.5)
  TEST_REAL_SIMILAR(fa.intensity, 678.0)

  TEST_EQUAL(toAnnotation_("y5+").charge, 1)
  TEST_STRING_EQUAL(toAnnotation_("y5+").annotation, "y5+")

  // negative mode, e.g. nucleic acid fragments
  TEST_EQUAL(toAnnotation_("c1-").charge, -1)
  TEST_STRING_EQUAL(toAnnotation_("c1-").annotation, "c1-")
  TEST_EQUAL(toAnnotation_("a3-B--").charge, -2)

  // a neutral loss before the charge, and the other two notations a name can use
  TEST_EQUAL(toAnnotation_("y3-H2O1++").charge, 2)
  TEST_STRING_EQUAL(toAnnotation_("y3-H2O1++").annotation, "y3-H2O1++")
  TEST_EQUAL(toAnnotation_("y5+3").charge, 3)
  TEST_EQUAL(toAnnotation_("b2^2").charge, 2)

  // a name that spells no charge leaves 'charge' at 0, which is distinguishable from charge 1
  TEST_EQUAL(toAnnotation_("y5").charge, 0)
  TEST_STRING_EQUAL(toAnnotation_("y5").annotation, "y5")
  // ... and a neutral loss ending in a digit is not a charge
  TEST_EQUAL(toAnnotation_("y1-H2O1").charge, 0)

  // surrounding whitespace is not part of the annotation
  TEST_STRING_EQUAL(toAnnotation_("  y5++  ").annotation, "y5++")
}
END_SECTION

START_SECTION(([EXTRA] toPeakAnnotation() never rewrites the ion name the producer wrote))
{
  // The identification view synchronises annotations back into the PeptideHit on every spectrum
  // switch, not only when the user edits something. So whatever the view does to show the charge must
  // not reach the name, or merely looking at a spectrum rewrites the stored data.
  // NucleicAcidSearchEngine writes "a2-B" with the charge in PeakAnnotation::charge alone; that name
  // has to survive untouched, otherwise it becomes "a2-B-" and the sequence diagram files the ion
  // under "a" instead of "a-B".
  PeptideHit::PeakAnnotation fa = toAnnotation_("a2-B", -1);
  TEST_STRING_EQUAL(fa.annotation, "a2-B")
  TEST_EQUAL(fa.charge, -1)

  // same for the cross-link names, which also leave the charge to the member
  TEST_STRING_EQUAL(toAnnotation_("[alpha|ci$y3]", 2).annotation, "[alpha|ci$y3]")
  TEST_EQUAL(toAnnotation_("[alpha|ci$y3]", 2).charge, 2)

  // and for an mzPAF name, which must stay parseable as mzPAF
  TEST_STRING_EQUAL(toAnnotation_("y3", 1).annotation, "y3")
  TEST_EQUAL(toAnnotation_("y3", 1).charge, 1)

  // a name that spells the charge keeps it and is still returned unchanged
  TEST_STRING_EQUAL(toAnnotation_("y3+", 1).annotation, "y3+")
  TEST_EQUAL(toAnnotation_("y3+", 1).charge, 1)

  // when the two disagree the name wins: it is what the producer actually wrote
  TEST_EQUAL(toAnnotation_("y3+", 2).charge, 1)

  // an unknown charge stays unknown rather than being invented
  TEST_EQUAL(toAnnotation_("y3", 0).charge, 0)
}
END_SECTION

START_SECTION(([EXTRA] toPeakAnnotation() keeps the free text below the ion name))
{
  // only the first line is the ion name, the rest is whatever the user typed and has to survive
  PeptideHit::PeakAnnotation fa = toAnnotation_("y5++\nmy comment");
  TEST_STRING_EQUAL(fa.annotation, "y5++\nmy comment")
  TEST_EQUAL(fa.charge, 2)

  // more than one extra line (older code kept only the first of them)
  TEST_STRING_EQUAL(toAnnotation_("y5++\nline one\nline two").annotation, "y5++\nline one\nline two")

  // a blank line the user put between the ion name and their comment is text, not a separator
  TEST_STRING_EQUAL(toAnnotation_("y5++\n\nmy comment").annotation, "y5++\n\nmy comment")
  TEST_EQUAL(toAnnotation_("y5++\n\nmy comment").charge, 2)

  // line endings are normalised to '\n', a CRLF is one line break and not an empty line
  TEST_STRING_EQUAL(toAnnotation_("y5++\r\nmy comment").annotation, "y5++\nmy comment")
  // a lone CR is a line break too, the same way the identification view treats it when it builds
  // the label -- it used to be deleted here, which joined the two lines into one
  TEST_STRING_EQUAL(toAnnotation_("y5++\rmy comment").annotation, "y5++\nmy comment")
  TEST_EQUAL(toAnnotation_("y5++\rmy comment").charge, 2)

  // the charge is read from the ion name, never from the free text below it
  TEST_EQUAL(toAnnotation_("y5\ncharge is 2+").charge, 0)
  TEST_STRING_EQUAL(toAnnotation_("y5\ncharge is 2+").annotation, "y5\ncharge is 2+")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
