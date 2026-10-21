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

#include <OpenMS/ANALYSIS/QUANTITATION/DDAWorkflowCommons.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>

#include <map>
#include <string>

///////////////////////////

START_TEST(DDAWorkflowCommons, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((static std::map<std::string, std::string> mapId2MzMLs(const std::map<std::string, std::string>& m2i)))
{
  // mapId2MzMLs simply inverts the mzML->id map into an id->mzML map
  std::map<std::string, std::string> m2i;
  m2i["/a/runA.mzML"] = "/x/runA.idXML";
  m2i["/b/runB.mzML"] = "/y/runB.idXML";

  std::map<std::string, std::string> i2m = DDAWorkflowCommons::mapId2MzMLs(m2i);

  TEST_EQUAL(i2m.size(), 2)
  TEST_STRING_EQUAL(i2m["/x/runA.idXML"], "/a/runA.mzML")
  TEST_STRING_EQUAL(i2m["/y/runB.idXML"], "/b/runB.mzML")

  // empty input -> empty output
  std::map<std::string, std::string> empty_in;
  TEST_EQUAL(DDAWorkflowCommons::mapId2MzMLs(empty_in).empty(), true)
}
END_SECTION

START_SECTION((static std::map<std::string, std::string> mapMzML2Ids(StringList& in, StringList& in_ids)))
{
  // Matching base names in the same order -> map of absolute mzML path -> absolute id path
  StringList in{"dirA/sample1.mzML", "dirB/sample2.mzML"};
  StringList ids{"idsA/sample1.idXML", "idsB/sample2.idXML"};

  std::map<std::string, std::string> mz2id = DDAWorkflowCommons::mapMzML2Ids(in, ids);

  TEST_EQUAL(mz2id.size(), 2)
  TEST_STRING_EQUAL(mz2id[File::absolutePath("dirA/sample1.mzML")], File::absolutePath("idsA/sample1.idXML"))
  TEST_STRING_EQUAL(mz2id[File::absolutePath("dirB/sample2.mzML")], File::absolutePath("idsB/sample2.idXML"))

  // Different number of files -> rejected
  StringList in_short{"a.mzML"};
  StringList ids_long{"a.idXML", "b.idXML"};
  TEST_EXCEPTION(Exception::IllegalArgument, DDAWorkflowCommons::mapMzML2Ids(in_short, ids_long))

  // Same base names but different order -> rejected
  StringList in_ord{"a.mzML", "b.mzML"};
  StringList ids_ord{"b.idXML", "a.idXML"};
  TEST_EXCEPTION(Exception::IllegalArgument, DDAWorkflowCommons::mapMzML2Ids(in_ord, ids_ord))
}
END_SECTION

START_SECTION(([EXTRA] FAIMS: calculateSeeds and estimateMedianChromatographicFWHM per compensation voltage))
{
  // Two peptides, each eluting at one of two interleaved FAIMS compensation voltages: MS1 spectra alternate between
  // -45 V and -60 V, one per second. Each peptide's isotope pattern (2+) has a Gaussian elution profile (sigma 6 s) in
  // the spectra of its CV only.
  const std::vector<std::pair<double, AASequence>> analytes = {
    {-45.0, AASequence::fromString("LVNELTEFAK")}, {-60.0, AASequence::fromString("YLYEIARR")}};
  const int charge = 2;
  const double apex = 60.0, sigma = 6.0;
  MSExperiment exp;
  for (int t = 0; t < 120; ++t)
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.setRT(static_cast<double>(t));
    const double cv = (t % 2 == 0) ? -45.0 : -60.0;
    spec.setDriftTime(cv);
    spec.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
    spec.setType(SpectrumSettings::SpectrumType::CENTROID);
    for (const auto& [analyte_cv, seq] : analytes)
    {
      if (analyte_cv != cv) continue;
      const double elution = 1e7 * std::exp(-0.5 * std::pow((t - apex) / sigma, 2));
      const IsotopeDistribution dist = CoarseIsotopePatternGenerator(4).estimateFromPeptideWeight(seq.getMonoWeight());
      for (Size i = 0; i < dist.size(); ++i)
      {
        spec.emplace_back(seq.getMZ(charge) + i * Constants::C13C12_MASSDIFF_U / charge, elution * dist[i].getIntensity());
      }
    }
    spec.sortByPosition();
    spec.setNativeID("scan=" + std::to_string(t + 1));
    exp.addSpectrum(std::move(spec));
  }
  exp.updateRanges();

  // The peak width is measured per CV: a trace of one CV is not broken up by the spectra of the other
  const double fwhm = DDAWorkflowCommons::estimateMedianChromatographicFWHM(exp);
  TEST_TRUE(fwhm > 0.0)

  FeatureMap seeds;
  DDAWorkflowCommons::calculateSeeds(exp, 0.0, seeds, 2.355 * sigma, "2:2", 3.0);
  TEST_TRUE(!seeds.empty())
  // every seed carries the CV it was found at, the one its peptide elutes at
  Size labelled = 0, wrong_cv = 0;
  bool found_minus_45 = false, found_minus_60 = false;
  for (const Feature& seed : seeds)
  {
    if (!seed.metaValueExists(Constants::UserParam::FAIMS_CV)) continue;
    ++labelled;
    const double cv = seed.getMetaValue(Constants::UserParam::FAIMS_CV);
    if (std::fabs(cv + 45.0) < 0.01) found_minus_45 = true;
    if (std::fabs(cv + 60.0) < 0.01) found_minus_60 = true;
    for (const auto& [analyte_cv, seq] : analytes)
    {
      if (std::fabs(seed.getMZ() - seq.getMZ(charge)) < 0.05 && cv != analyte_cv) ++wrong_cv;
    }
  }
  TEST_EQUAL(labelled, seeds.size())
  TEST_EQUAL(wrong_cv, 0)
  // both CV groups produced seeds
  TEST_TRUE(found_minus_45)
  TEST_TRUE(found_minus_60)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
