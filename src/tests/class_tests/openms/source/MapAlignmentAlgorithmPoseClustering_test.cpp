// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/ConversionHelper.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

using namespace std;
using namespace OpenMS;

// exposes the protected peak and feature selection
class TestAligner : public MapAlignmentAlgorithmPoseClustering
{
public:
  using MapAlignmentAlgorithmPoseClustering::selectPeaks_;
  using MapAlignmentAlgorithmPoseClustering::selectSpreadPeaks_;
  using MapAlignmentAlgorithmPoseClustering::selectFeatures_;
};

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmPoseClustering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmPoseClustering* ptr = nullptr;
MapAlignmentAlgorithmPoseClustering* nullPointer = nullptr;
START_SECTION((MapAlignmentAlgorithmPoseClustering()))
	ptr = new MapAlignmentAlgorithmPoseClustering();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithmPoseClustering()))
	delete ptr;
END_SECTION

START_SECTION((void setReference(const FeatureMap& map)))
{
  NOT_TESTABLE // tested together with "align"
}
END_SECTION

START_SECTION((void setReference(const PeakMap& map)))
{
  NOT_TESTABLE // tested together with "align"
}
END_SECTION

START_SECTION((void align(const PeakMap& map, TransformationDescription& trafo)))
{
  MzMLFile f;
  std::vector<PeakMap > maps(2);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in1.mzML.gz"), maps[0]);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in2.mzML.gz"), maps[1]);

  MapAlignmentAlgorithmPoseClustering aligner;
  aligner.setReference(maps[0]);

  TransformationDescription trafo;
  aligner.align(maps[1], trafo);

  TEST_EQUAL(trafo.getModelType(), "linear");
  TEST_EQUAL(trafo.getDataPoints().size(), 307);
  
  // @TODO: can we get the slope/intercept without fitting a model again?
  TransformationModelLinear lm(trafo.getDataPoints(),
                               trafo.getModelParameters());
  double slope, intercept;
  std::string x_weight, y_weight;
  double x_datum_min, x_datum_max, y_datum_min, y_datum_max;
  lm.getParameters(slope, intercept, x_weight, y_weight, x_datum_min, x_datum_max, y_datum_min, y_datum_max);
  TEST_REAL_SIMILAR(slope, 1.01164);
  TEST_REAL_SIMILAR(intercept, -32.0912);
}
END_SECTION

/////
// selectPeaks_() picks the peaks that setReference(PeakMap) and align(PeakMap) align

START_SECTION((static void selectPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n)))
{
  PeakMap mse;
  MSSpectrum mss;
  Peak1D p;
  for (UInt m = 0; m < 3; ++m)
  {
    mss.clear(true);
    for (UInt i = 0; i < 4; ++i)
    {
      p.setMZ(10 * m + i + 100.35);
      p.setIntensity(900 + 7 * m + 5 * i);
      mss.push_back(p);
    }
    mse.addSpectrum(mss);
    mse.getSpectra().back().setRT(m * 5);
  }

  ConsensusMap cm;
  TestAligner::selectPeaks_(33, mse, cm, 8);
  TEST_EQUAL(cm.size(), 8);
  TEST_EQUAL(cm.getColumnHeaders()[33].size, 8);
  TEST_EQUAL(cm.front().getIntensity(), 929);
  TEST_EQUAL(cm.back().getIntensity(), 912);
  for (const ConsensusFeature& cf : cm)
  {
    TEST_EQUAL(cf.size(), 1);
    TEST_EQUAL(cf.begin()->getMapIndex(), 33);
  }
}
END_SECTION

/////
// Inputs whose getSize() exceeds their MS1 peak count. get2DData() collects MS1 peaks
// only, whereas getSize() also counts the peaks of other MS levels and all chromatogram
// points, so n must be capped by the MS1 peak count. All MS2 and chromatogram intensities
// lie far above every MS1 intensity: anything leaking into the result would sort to the
// front and be caught by the content checks below.

// adds three MS2 spectra with five peaks each (15 peaks, intensities >= 1000)
auto add_ms2_spectra = [](PeakMap& map)
{
  MSSpectrum ms2;
  Peak1D p;
  for (UInt m = 0; m < 3; ++m)
  {
    ms2.clear(true);
    ms2.setMSLevel(2);
    ms2.setRT(10.0 * m + 1.0);
    for (UInt i = 0; i < 5; ++i)
    {
      p.setMZ(300.0 + i);
      p.setIntensity(1000.0f + 10 * m + i);
      ms2.push_back(p);
    }
    map.addSpectrum(ms2);
  }
};

// adds one chromatogram with ten points (intensities >= 5000)
auto add_chromatogram = [](PeakMap& map)
{
  MSChromatogram chrom;
  ChromatogramPeak cp;
  for (UInt i = 0; i < 10; ++i)
  {
    cp.setRT(1.0 * i);
    cp.setIntensity(5000.0f + i);
    chrom.push_back(cp);
  }
  map.addChromatogram(chrom);
};

// mixed input: two MS1 spectra with three peaks each (unique intensities 100, 101, 102 at
// RT 0 and 110, 111, 112 at RT 10), plus the MS2 spectra and the chromatogram from above.
// ms1_by_intensity holds the six MS1 peaks in the order selectPeaks_() must emit them.
PeakMap mixed;
std::vector<Peak2D> ms1_by_intensity;
{
  MSSpectrum ms1;
  Peak1D p;
  Peak2D p2;
  for (UInt m = 0; m < 2; ++m)
  {
    ms1.clear(true);
    ms1.setMSLevel(1);
    ms1.setRT(10.0 * m);
    for (UInt i = 0; i < 3; ++i)
    {
      p.setMZ(200.0 + 10 * m + i);
      p.setIntensity(100.0f + 10 * m + i);
      ms1.push_back(p);
      p2.setRT(ms1.getRT());
      p2.setMZ(p.getMZ());
      p2.setIntensity(p.getIntensity());
      ms1_by_intensity.push_back(p2);
    }
    mixed.addSpectrum(ms1);
  }
  add_ms2_spectra(mixed);
  add_chromatogram(mixed);
  std::sort(ms1_by_intensity.begin(), ms1_by_intensity.end(),
            [](const Peak2D& a, const Peak2D& b) { return a.getIntensity() > b.getIntensity(); });
}
const Size n_ms1 = ms1_by_intensity.size();
const UInt64 map_index = 5;

// checks that out holds exactly the 'expected' most intense MS1 peaks, most intense first,
// and that the column header reports the number of features written
auto check_top_ms1 = [&](ConsensusMap& out, Size expected)
{
  TEST_EQUAL(out.size(), expected)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, expected)
  for (Size i = 0; i < std::min(out.size(), expected); ++i)
  {
    TEST_REAL_SIMILAR(out[i].getIntensity(), ms1_by_intensity[i].getIntensity())
    TEST_REAL_SIMILAR(out[i].getRT(), ms1_by_intensity[i].getRT())
    TEST_REAL_SIMILAR(out[i].getMZ(), ms1_by_intensity[i].getMZ())
    TEST_EQUAL(out[i].size(), 1)
    TEST_EQUAL(out[i].begin()->getMapIndex(), map_index)
  }
};

START_SECTION([EXTRA] selectPeaks_ caps n by the number of MS1 peaks rather than by getSize())
{
  // the input really is in the problematic regime: far more peaks in total than MS1 peaks
  mixed.updateRanges();
  TEST_EQUAL(n_ms1, 6)
  TEST_EQUAL(mixed.getSize(), 6 + 15 + 10)

  ConsensusMap out;

  // n = Size(-1) ('max_num_peaks_considered' -1): all MS1 peaks and nothing else
  TestAligner::selectPeaks_(map_index, mixed, out, Size(-1));
  check_top_ms1(out, n_ms1);

  // n above the MS1 peak count but below getSize()
  TestAligner::selectPeaks_(map_index, mixed, out, 20);
  check_top_ms1(out, n_ms1);

  // n below the MS1 peak count still selects the n most intense MS1 peaks
  TestAligner::selectPeaks_(map_index, mixed, out, 4);
  check_top_ms1(out, 4);
}
END_SECTION

START_SECTION([EXTRA] selectPeaks_ without MS1 peaks yields an empty ConsensusMap)
{
  // MS2-only input; n = number of spectra and n = Size(-1) (all peaks)
  PeakMap ms2_only;
  add_ms2_spectra(ms2_only);
  TEST_EQUAL(ms2_only.getNrSpectra(), 3)
  TEST_EQUAL(ms2_only.getSize(), 15)

  ConsensusMap out;
  TestAligner::selectPeaks_(map_index, ms2_only, out, ms2_only.size());
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)

  TestAligner::selectPeaks_(map_index, ms2_only, out, Size(-1));
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)

  // chromatogram points are not peaks either
  PeakMap chrom_only;
  add_chromatogram(chrom_only);
  TEST_EQUAL(chrom_only.getSize(), 10)

  TestAligner::selectPeaks_(map_index, chrom_only, out, Size(-1));
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)
}
END_SECTION

START_SECTION([EXTRA] selectPeaks_ orders equal intensities by RT and m/z)
{
  // many peaks with the same intensity: which of them are kept (and in which order) must not depend on
  // the standard library's partial_sort
  PeakMap ties;
  for (int s = 0; s < 3; ++s)
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.setRT(30.0 - 10.0 * s); // RT 30, 20, 10
    for (int p = 0; p < 20; ++p)
    {
      spec.push_back(Peak1D(900.0 - 10.0 * p, 5.0f)); // m/z 900, 890, ..., 710
    }
    ties.addSpectrum(spec);
  }
  ties[1][7].setIntensity(9.0f); // RT 20, m/z 830 is the most intense peak

  ConsensusMap out;
  TestAligner::selectPeaks_(map_index, ties, out, 4);
  TEST_EQUAL(out.size(), 4)
  ABORT_IF(out.size() != 4)
  TEST_REAL_SIMILAR(out[0].getRT(), 20.0)
  TEST_REAL_SIMILAR(out[0].getMZ(), 830.0)
  TEST_REAL_SIMILAR(out[1].getRT(), 10.0)
  TEST_REAL_SIMILAR(out[1].getMZ(), 710.0)
  TEST_REAL_SIMILAR(out[2].getRT(), 10.0)
  TEST_REAL_SIMILAR(out[2].getMZ(), 720.0)
  TEST_REAL_SIMILAR(out[3].getRT(), 10.0)
  TEST_REAL_SIMILAR(out[3].getMZ(), 730.0)
}
END_SECTION

START_SECTION([EXTRA] selectPeaks_ rejects NaN intensities)
{
  PeakMap nan_map;
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.push_back(Peak1D(500.0, 10.0f));
  spec.push_back(Peak1D(600.0, std::numeric_limits<float>::quiet_NaN()));
  nan_map.addSpectrum(spec);
  ConsensusMap out;
  TEST_EXCEPTION(Exception::InvalidValue, TestAligner::selectPeaks_(map_index, nan_map, out, 1))
  TEST_EXCEPTION(Exception::InvalidValue, TestAligner::selectPeaks_(map_index, nan_map, out, Size(-1)))

  // NaN outside MS1 does not matter
  nan_map[0].setMSLevel(2);
  TestAligner::selectPeaks_(map_index, nan_map, out, Size(-1));
  TEST_EQUAL(out.size(), 0)
}
END_SECTION

/////
// selectSpreadPeaks_() picks one peak per analyte, spread over the map

START_SECTION((static void selectSpreadPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n, double rt_window, double mz_tolerance_ppm, Size max_charge, Size rt_bins, Size mz_bins)))
{
  // analyte A at m/z 500 elutes over RT 0..9 (apex 1000 at RT 5) with its first isotope (60%); analyte B at m/z 800
  // elutes at the same time (apex 300); analyte C has A's m/z but elutes much later (RT 100)
  const double isotope = Constants::C13C12_MASSDIFF_U;
  PeakMap elution;
  for (int s = 0; s < 10; ++s)
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.setRT(s);
    const float a = 1000.0f - 100.0f * std::abs(s - 5);
    spec.push_back(Peak1D(500.0, a));
    spec.push_back(Peak1D(500.0 + isotope, 0.6f * a));
    spec.push_back(Peak1D(800.0, 300.0f - 30.0f * std::abs(s - 5)));
    elution.addSpectrum(spec);
  }
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.setRT(100.0);
    spec.push_back(Peak1D(500.0, 200.0f));
    elution.addSpectrum(spec);
  }

  // the most intense peaks are all from A
  ConsensusMap out;
  TestAligner::selectPeaks_(map_index, elution, out, 3);
  TEST_EQUAL(out.size(), 3)
  ABORT_IF(out.size() != 3)
  TEST_REAL_SIMILAR(out[0].getMZ(), 500.0)
  TEST_REAL_SIMILAR(out[1].getMZ(), 500.0)
  TEST_REAL_SIMILAR(out[2].getMZ(), 500.0)

  // one peak (the apex) per analyte
  TestAligner::selectSpreadPeaks_(map_index, elution, out, 3, 20.0, 10.0, 2, 1, 1);
  TEST_EQUAL(out.size(), 3)
  ABORT_IF(out.size() != 3)
  TEST_REAL_SIMILAR(out[0].getRT(), 5.0)
  TEST_REAL_SIMILAR(out[0].getMZ(), 500.0)
  TEST_REAL_SIMILAR(out[0].getIntensity(), 1000.0)
  TEST_REAL_SIMILAR(out[1].getRT(), 5.0)
  TEST_REAL_SIMILAR(out[1].getMZ(), 800.0)
  TEST_REAL_SIMILAR(out[2].getRT(), 100.0)
  TEST_REAL_SIMILAR(out[2].getMZ(), 500.0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 3)
  TEST_EQUAL(out[0].size(), 1)
  TEST_EQUAL(out[0].begin()->getMapIndex(), map_index)

  // fewer analytes than n: fewer peaks
  TestAligner::selectSpreadPeaks_(map_index, elution, out, 10, 20.0, 10.0, 2, 1, 1);
  TEST_EQUAL(out.size(), 3)

  // without isotope suppression, A's isotope counts as an analyte of its own
  TestAligner::selectSpreadPeaks_(map_index, elution, out, 10, 20.0, 10.0, 0, 1, 1);
  TEST_EQUAL(out.size(), 4)
  ABORT_IF(out.size() != 4)
  TEST_REAL_SIMILAR(out[1].getMZ(), 500.0 + isotope)

  // isotopes of charge 2 are half as far apart: suppressed with max_charge 2, not with max_charge 1
  PeakMap charge2;
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.push_back(Peak1D(600.0, 1000.0f));
    spec.push_back(Peak1D(600.0 + isotope / 2.0, 800.0f));
    charge2.addSpectrum(spec);
  }
  TestAligner::selectSpreadPeaks_(map_index, charge2, out, 5, 20.0, 10.0, 2, 1, 1);
  TEST_EQUAL(out.size(), 1)
  TestAligner::selectSpreadPeaks_(map_index, charge2, out, 5, 20.0, 10.0, 1, 1, 1);
  TEST_EQUAL(out.size(), 2)

  // spreading: four intense analytes early (RT 10..40), two weak ones late (RT 70, 90)
  PeakMap early_late;
  const std::vector<std::pair<double, std::pair<double, float>>> analytes{
    {10.0, {300.0, 1000.0f}}, {20.0, {400.0, 900.0f}}, {30.0, {500.0, 800.0f}}, {40.0, {600.0, 700.0f}},
    {70.0, {350.0, 20.0f}}, {90.0, {450.0, 10.0f}}};
  for (const auto& analyte : analytes)
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.setRT(analyte.first);
    spec.push_back(Peak1D(analyte.second.first, analyte.second.second));
    early_late.addSpectrum(spec);
  }
  // one RT bin: the four most intense (all early)
  TestAligner::selectSpreadPeaks_(map_index, early_late, out, 4, 5.0, 10.0, 2, 1, 1);
  TEST_EQUAL(out.size(), 4)
  ABORT_IF(out.size() != 4)
  TEST_REAL_SIMILAR(out[3].getIntensity(), 700.0)
  // two RT bins: two from each half, round-robin
  TestAligner::selectSpreadPeaks_(map_index, early_late, out, 4, 5.0, 10.0, 2, 2, 1);
  TEST_EQUAL(out.size(), 4)
  ABORT_IF(out.size() != 4)
  TEST_REAL_SIMILAR(out[0].getIntensity(), 1000.0)
  TEST_REAL_SIMILAR(out[1].getIntensity(), 900.0)
  TEST_REAL_SIMILAR(out[2].getIntensity(), 20.0)
  TEST_REAL_SIMILAR(out[3].getIntensity(), 10.0)

  // empty, MS2 only, n = 0
  TestAligner::selectSpreadPeaks_(map_index, PeakMap(), out, 5, 20.0, 10.0, 2, 10, 10);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)
  PeakMap ms2;
  {
    MSSpectrum spec;
    spec.setMSLevel(2);
    spec.push_back(Peak1D(500.0, 5.0f));
    ms2.addSpectrum(spec);
  }
  TestAligner::selectSpreadPeaks_(map_index, ms2, out, 5, 20.0, 10.0, 2, 10, 10);
  TEST_EQUAL(out.size(), 0)
  TestAligner::selectSpreadPeaks_(map_index, early_late, out, 0, 20.0, 10.0, 2, 10, 10);
  TEST_EQUAL(out.size(), 0)

  // NaN intensities are rejected
  PeakMap nan_map;
  {
    MSSpectrum spec;
    spec.setMSLevel(1);
    spec.push_back(Peak1D(500.0, std::numeric_limits<float>::quiet_NaN()));
    nan_map.addSpectrum(spec);
  }
  TEST_EXCEPTION(Exception::InvalidValue, TestAligner::selectSpreadPeaks_(map_index, nan_map, out, 5, 20.0, 10.0, 2, 1, 1))
}
END_SECTION

START_SECTION([EXTRA] align(PeakMap) with peak_selection:method 'spread')
{
  MzMLFile f;
  std::vector<PeakMap> maps(2);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in1.mzML.gz"), maps[0]);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in2.mzML.gz"), maps[1]);

  MapAlignmentAlgorithmPoseClustering aligner;
  Param param = aligner.getParameters();
  param.setValue("peak_selection:method", "spread");
  aligner.setParameters(param);
  aligner.setReference(maps[0]);

  TransformationDescription trafo;
  aligner.align(maps[1], trafo);
  TEST_EQUAL(trafo.getModelType(), "linear");
  STATUS("data points with 'spread': " << trafo.getDataPoints().size());
  TEST_EQUAL(trafo.getDataPoints().empty(), false);
}
END_SECTION

/////
// selectFeatures_() picks the features that setReference(FeatureMap) and align(FeatureMap) align

START_SECTION((static void selectFeatures_(UInt64 map_index, const FeatureMap& map, ConsensusMap& output, Size n)))
{
  // sorted by RT, the most intense features come last
  FeatureMap fm;
  fm.setUniqueId(42);
  for (UInt i = 0; i < 5; ++i)
  {
    Feature f;
    f.setRT(10.0 * i);
    f.setMZ(400.0 + i);
    f.setIntensity(100.0f * (i + 1));
    f.setUniqueId(i + 1);
    fm.push_back(f);
  }
  fm.getProteinIdentifications().resize(1);

  ConsensusMap out;
  // the two most intense features (not the first two), in map order
  TestAligner::selectFeatures_(map_index, fm, out, 2);
  TEST_EQUAL(out.size(), 2)
  ABORT_IF(out.size() != 2)
  TEST_REAL_SIMILAR(out[0].getRT(), 30.0)
  TEST_REAL_SIMILAR(out[1].getRT(), 40.0)
  TEST_EQUAL(out[0].begin()->getUniqueId(), 4)
  TEST_EQUAL(out[1].begin()->getUniqueId(), 5)
  TEST_EQUAL(out[0].begin()->getMapIndex(), map_index)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 5) // the input size, as in MapConversion::convert()
  TEST_EQUAL(out.getUniqueId(), 42)
  TEST_EQUAL(out.getProteinIdentifications().size(), 1)

  // equal intensities: earlier RT first
  fm[3].setIntensity(500.0f); // ties with fm[4]
  TestAligner::selectFeatures_(map_index, fm, out, 1);
  TEST_EQUAL(out.size(), 1)
  ABORT_IF(out.size() != 1)
  TEST_REAL_SIMILAR(out[0].getRT(), 30.0)

  // all features: the same as MapConversion::convert()
  ConsensusMap converted;
  MapConversion::convert(map_index, fm, converted);
  TestAligner::selectFeatures_(map_index, fm, out, Size(-1));
  TEST_EQUAL(out.size(), converted.size())
  for (Size i = 0; i < std::min(out.size(), converted.size()); ++i)
  {
    TEST_EQUAL(out[i], converted[i])
  }

  // NaN intensities are rejected
  fm[2].setIntensity(std::numeric_limits<float>::quiet_NaN());
  TEST_EXCEPTION(Exception::InvalidValue, TestAligner::selectFeatures_(map_index, fm, out, 2))
  TEST_EXCEPTION(Exception::InvalidValue, TestAligner::selectFeatures_(map_index, fm, out, Size(-1)))
}
END_SECTION

START_SECTION((void align(const FeatureMap& map, TransformationDescription& trafo)))
{
  // Tested extensively in TEST/TOPP
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((void align(const ConsensusMap& map, TransformationDescription& trafo)))
{
  // Tested extensively in TEST/TOPP
  NOT_TESTABLE;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
