// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/CheckedCast.h>
#include <OpenMS/DATASTRUCTURES/RegularExpression.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <vector>

int main()
{
  OpenMS::MetaInfo meta;
  meta.setValue("answer", 42);
  if (meta.getValue("answer") != OpenMS::DataValue(42)) return 1;
  OpenMS::IdentificationData ids;
  auto file = ids.registerInputFile(OpenMS::IdentificationData::InputFile("test"));
  ids.registerObservation(OpenMS::IdentificationData::Observation("scan=42", file));
  if (ids.getObservations().size() != 1) return 2;
  const OpenMS::RegularExpression regex(R"(scan=(?<SCAN>\d+))");
  if (OpenMS::SpectrumNativeIDParser::extractScanNumber("scan=42", regex) != 42) return 3;
  auto copy = regex;
  copy.assign("other");
  if (! regex.search("scan=42") || copy.search("scan=42")) return 4;
  OpenMS::Math::RandomShuffler shuffler(42);
  std::vector<int> values {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  shuffler.portable_random_shuffle(values.begin(), values.end());
  if (values != std::vector<int> {2, 6, 5, 4, 8, 3, 0, 11, 1, 10, 7, 9}) return 5;
  if (OpenMS::checkedCast<int>(values.size()) != 12) return 6;
  try
  {
    (void)OpenMS::checkedCast<unsigned>(-1);
    return 7;
  }
  catch (const OpenMS::Exception::OutOfRange&)
  {
  }
  return OpenMS::Math::binomial_cdf_complement(2, 1, 0.5) == 0.75 ? 0 : 8;
}
