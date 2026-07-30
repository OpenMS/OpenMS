// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>

#include <optional>

using namespace OpenMS;

START_TEST(Annotation1DPeakItem, "$Id$")

START_SECTION((PeptideHit::PeakAnnotation toPeakAnnotation() const))
{
  Annotation1DPeakItem<Peak1D> item(
    Peak1D(499.999, 1000.0f), "y4+", QColor(), std::optional<double>(500.0));

  PeptideHit::PeakAnnotation annotation = item.toPeakAnnotation();
  TEST_EQUAL(annotation.annotation, "y4")
  TEST_EQUAL(annotation.charge, 1)
  TEST_REAL_SIMILAR(annotation.mz, 499.999)
  TEST_REAL_SIMILAR(annotation.intensity, 1000.0)
  TEST_TRUE(annotation.theoretical_mz.has_value())
  TEST_REAL_SIMILAR(*annotation.theoretical_mz, 500.0)

  item.clearTheoreticalMZ();
  annotation = item.toPeakAnnotation();
  TEST_FALSE(annotation.theoretical_mz.has_value())
}
END_SECTION

END_TEST
