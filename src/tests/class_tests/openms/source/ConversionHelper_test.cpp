// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConversionHelper.h>
///////////////////////////

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((template < typename FeatureT > static void convert(UInt64 const input_map_index, FeatureMap< FeatureT > const &input_map, ConsensusMap &output_map, Size n=-1)))
{

  FeatureMap fm;
  Feature f;
  for ( UInt i = 0; i < 3; ++i )
  {
    f.setRT(i*77.7);
    f.setMZ(i+100.35);
    f.setUniqueId(i*33+17);
    fm.push_back(f);
  }
  ConsensusMap cm;
  MapConversion::convert(33,fm,cm);

  TEST_EQUAL(cm.size(),3);
  TEST_EQUAL(cm.getColumnHeaders()[33].size,3);
  for ( UInt i = 0; i < 3; ++i )
  {
    TEST_EQUAL(cm[i].size(),1);
    TEST_EQUAL(cm[i].begin()->getMapIndex(),33);
    TEST_EQUAL(cm[i].begin()->getUniqueId(),i*33+17);
    TEST_REAL_SIMILAR(cm[i].begin()->getRT(),i*77.7);
    TEST_REAL_SIMILAR(cm[i].begin()->getMZ(),i+100.35);
  }

cm.clear();
MapConversion::convert(33,fm,cm,2);
TEST_EQUAL(cm.size(),2);
TEST_EQUAL(cm.getColumnHeaders()[33].size,3);

}
END_SECTION

/////

// Prepare data
PeakMap mse;
{
  MSSpectrum mss;
  Peak1D p;
  for ( UInt m = 0; m < 3; ++m )
  {
    mss.clear(true);
    for ( UInt i = 0; i < 4; ++i )
    {
      p.setMZ( 10* m + i + 100.35);
      p.setIntensity( 900 + 7*m + 5*i );
      mss.push_back(p);
    }
    mse.addSpectrum(mss);
    mse.getSpectra().back().setRT(m*5);
  }
}

START_SECTION((static void convert(UInt64 const input_map_index, PeakMap & input_map, ConsensusMap& output_map, Size n = -1)))
{

  ConsensusMap cm;

  MapConversion::convert(33,mse,cm,8);

  TEST_EQUAL(cm.size(),8);

  for ( UInt i = 0; i < cm.size(); ++i)
  {
    STATUS("\n" << i << ": " << cm[i] );
  }

  TEST_EQUAL(cm.back().getIntensity(),912);

}
END_SECTION

/////

ConsensusMap cm;
MapConversion::convert(33,mse,cm,8);

START_SECTION((template < typename FeatureT > static void convert(ConsensusMap const &input_map, const bool keep_uids, FeatureMap< FeatureT > &output_map)))
{
    FeatureMap out_fm;
    MapConversion::convert(cm, true, out_fm);

    TEST_EQUAL(cm.getUniqueId(), out_fm.getUniqueId());
    TEST_EQUAL(cm.getProteinIdentifications().size(), out_fm.getProteinIdentifications().size());
    TEST_EQUAL(cm.getUnassignedPeptideIdentifications().size(), out_fm.getUnassignedPeptideIdentifications().size());
    TEST_EQUAL(cm.size(), out_fm.size());

    for (Size i = 0; i < cm.size(); ++i)
    {
        TEST_EQUAL(cm[i], out_fm[i]);
    }

    out_fm.clear();
    MapConversion::convert(cm, false, out_fm);
    TEST_NOT_EQUAL(cm.getUniqueId(), out_fm.getUniqueId());

    for (Size i = 0; i < cm.size(); ++i)
    {
        TEST_REAL_SIMILAR(cm[i].getRT(), out_fm[i].getRT());
        TEST_REAL_SIMILAR(cm[i].getMZ(), out_fm[i].getMZ());
        TEST_REAL_SIMILAR(cm[i].getIntensity(), out_fm[i].getIntensity());

        TEST_NOT_EQUAL(cm[i].getUniqueId(), out_fm[i].getUniqueId());
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
// ms1_by_intensity holds the six MS1 peaks in the order convert() must emit them.
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

START_SECTION([EXTRA] convert(PeakMap) caps n by the number of MS1 peaks, not by getSize())
{
  // the input really is in the problematic regime: far more peaks in total than MS1 peaks
  mixed.updateRanges();
  TEST_EQUAL(n_ms1, 6)
  TEST_EQUAL(mixed.getSize(), 6 + 15 + 10)

  ConsensusMap out;

  // default n (Size(-1)): all MS1 peaks and nothing else
  MapConversion::convert(map_index, mixed, out);
  check_top_ms1(out, n_ms1);

  // n above the MS1 peak count but below getSize()
  MapConversion::convert(map_index, mixed, out, 20);
  check_top_ms1(out, n_ms1);

  // n below the MS1 peak count still selects the n most intense MS1 peaks
  MapConversion::convert(map_index, mixed, out, 4);
  check_top_ms1(out, 4);
}
END_SECTION

START_SECTION([EXTRA] convert(PeakMap) without MS1 peaks yields an empty ConsensusMap)
{
  // MS2-only input; n = number of spectra (as FileConverter passes it) and the default n
  PeakMap ms2_only;
  add_ms2_spectra(ms2_only);
  TEST_EQUAL(ms2_only.getNrSpectra(), 3)
  TEST_EQUAL(ms2_only.getSize(), 15)

  ConsensusMap out;
  MapConversion::convert(map_index, ms2_only, out, ms2_only.size());
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)

  MapConversion::convert(map_index, ms2_only, out);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)

  // chromatogram points are not peaks either
  PeakMap chrom_only;
  add_chromatogram(chrom_only);
  TEST_EQUAL(chrom_only.getSize(), 10)

  MapConversion::convert(map_index, chrom_only, out);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getColumnHeaders()[map_index].size, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



