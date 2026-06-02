// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakMapExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;
using namespace std;

START_TEST(PeakMapExtractor, "$Id$")

START_SECTION(void extractPeakMaps(const OpenSwath::SpectrumAccessPtr&, std::vector<ExtractedPeakMap>&, const std::vector<ExtractionCoordinates>&, double, bool, double, const String&))
{
  typedef OpenMS::DataArrays::FloatDataArray FloatDataArray;

  std::shared_ptr<PeakMap> exp(new PeakMap);
  for (int scan = 0; scan < 2; ++scan)
  {
    MSSpectrum spectrum;
    spectrum.setRT(scan == 0 ? 10.0 : 20.0);

    FloatDataArray im_array;
    im_array.setName(Constants::UserParam::ION_MOBILITY);

    const double intensities[] = {10.0, 11.0, 12.0, 13.0, 14.0};
    const double mzs[] = {499.99, 500.00, 500.02, 500.03, 500.10};
    const double ims[] = {1.00, 1.05, 1.10, 1.15, 1.30};
    for (Size i = 0; i < 5; ++i)
    {
      Peak1D peak;
      peak.setMZ(mzs[i]);
      peak.setIntensity(intensities[i] + 100.0 * scan);
      spectrum.push_back(peak);
      im_array.push_back(ims[i]);
    }
    spectrum.getFloatDataArrays().push_back(im_array);
    exp->addSpectrum(spectrum);
  }

  OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  PeakMapExtractor extractor;
  std::vector<PeakMapExtractor::ExtractionCoordinates> coordinates;
  PeakMapExtractor::ExtractionCoordinates coord;
  coord.mz = 500.02;
  coord.ion_mobility = 1.10;
  coord.rt_start = 0.0;
  coord.rt_end = -1.0;
  coord.id = "tr1";
  coordinates.push_back(coord);

  coord.id = "tr2";
  coord.rt_start = 15.0;
  coord.rt_end = 25.0;
  coordinates.push_back(coord);

  std::vector<PeakMapExtractor::ExtractedPeakMap> output;
  extractor.extractPeakMaps(expptr, output, coordinates, 0.05, false, 0.20, "tophat");

  TEST_EQUAL(output.size(), 2)
  TEST_STRING_EQUAL(output[0].native_id, "tr1")
  TEST_EQUAL(output[0].mz.size(), 6)
  TEST_EQUAL(output[0].mz.size(), output[0].rt.size())
  TEST_EQUAL(output[0].mz.size(), output[0].ion_mobility.size())
  TEST_EQUAL(output[0].mz.size(), output[0].intensity.size())
  TEST_REAL_SIMILAR(output[0].rt[0], 10.0)
  TEST_REAL_SIMILAR(output[0].rt[3], 20.0)
  TEST_REAL_SIMILAR(output[0].mz[0], 500.00)
  TEST_REAL_SIMILAR(output[0].mz[2], 500.03)
  TEST_REAL_SIMILAR(output[0].intensity[0], 11.0)
  TEST_REAL_SIMILAR(output[0].intensity[5], 113.0)

  TEST_STRING_EQUAL(output[1].native_id, "tr2")
  TEST_EQUAL(output[1].mz.size(), 3)
  TEST_REAL_SIMILAR(output[1].rt[0], 20.0)
  TEST_REAL_SIMILAR(output[1].intensity[0], 111.0)
}
END_SECTION

START_SECTION([EXTRA] void extractPeakMaps throws without ion mobility data)
{
  std::shared_ptr<PeakMap> exp(new PeakMap);
  MSSpectrum spectrum;
  spectrum.setRT(10.0);
  spectrum.push_back(Peak1D(500.0, 100.0));
  exp->addSpectrum(spectrum);

  OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  PeakMapExtractor extractor;
  std::vector<PeakMapExtractor::ExtractionCoordinates> coordinates(1);
  coordinates[0].mz = 500.0;
  coordinates[0].ion_mobility = 1.0;
  coordinates[0].id = "tr1";

  std::vector<PeakMapExtractor::ExtractedPeakMap> output;
  TEST_EXCEPTION(Exception::IllegalArgument,
                 extractor.extractPeakMaps(expptr, output, coordinates, 0.05, false, 0.20, "tophat"))
}
END_SECTION

END_TEST
