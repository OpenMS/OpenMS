// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Holger Franken, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MassTrace.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

PeakType fillPeak(double rt, double mz, double it)
{
  PeakType p;
  p.setIntensity((OpenMS::Peak2D::IntensityType)it);
  p.setMZ(mz);
  p.setRT((rt));
  return p;
}

START_TEST(MassTrace, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassTrace* d_ptr = nullptr;
MassTrace* nullPointer = nullptr;
START_SECTION((MassTrace()))
{
    d_ptr = new MassTrace;
    TEST_NOT_EQUAL(d_ptr, nullPointer);
}
END_SECTION

START_SECTION((~MassTrace()))
{
    delete d_ptr;
}
END_SECTION


std::vector<PeakType> peak_vec;
std::list<PeakType> peak_lst;

PeakType tmp_peak0 = fillPeak(152.22, 230.10223, 542.0);
peak_vec.push_back(tmp_peak0);
peak_lst.push_back(tmp_peak0);

PeakType tmp_peak1 = fillPeak(153.23, 230.10235, 542293.0);
peak_vec.push_back(tmp_peak1);
peak_lst.push_back(tmp_peak1);

PeakType tmp_peak2 = fillPeak(154.21, 230.10181, 18282393.0);
peak_vec.push_back(tmp_peak2);
peak_lst.push_back(tmp_peak2);

PeakType tmp_peak3 = fillPeak(155.24, 230.10229, 33329535.0);
peak_vec.push_back(tmp_peak3);
peak_lst.push_back(tmp_peak3);

PeakType tmp_peak4 = fillPeak(156.233, 230.10116, 17342933.0);
peak_vec.push_back(tmp_peak4);
peak_lst.push_back(tmp_peak4);

PeakType tmp_peak5 = fillPeak(157.24, 230.10198, 333291.0);
peak_vec.push_back(tmp_peak5);
peak_lst.push_back(tmp_peak5);

PeakType tmp_peak6 = fillPeak(158.238, 230.10254, 339.0);
peak_vec.push_back(tmp_peak6);
peak_lst.push_back(tmp_peak6);

String si, sm, sr;
for (Size i = 0; i<peak_vec.size(); ++i)
{
  si += ", " + String(peak_vec[i].getIntensity());
  sm += ", " + String(peak_vec[i].getMZ());
  sr += ", " + String(peak_vec[i].getRT());
}
std::cout << sr << "\n" << sm << "\n" << si << "\n\n";


/////////////////////////////////////////////////////////////
// detailed constructors test
/////////////////////////////////////////////////////////////

START_SECTION((MassTrace(const std::list<PeakType>& trace_peaks)))
{
  MassTrace tmp_mt(peak_lst);

  std::list<PeakType>::const_iterator l_it = peak_lst.begin();

  for (MassTrace::const_iterator m_it = tmp_mt.begin(); m_it != tmp_mt.end(); ++m_it)
  {
      TEST_EQUAL(*l_it, *m_it);
      ++l_it;
  }

  TEST_REAL_SIMILAR(tmp_mt.getAverageMS1CycleTime(), (tmp_peak6.getRT() - tmp_peak0.getRT()) / (7-1));
}
END_SECTION

/////

START_SECTION((MassTrace(const std::vector<PeakType>& trace_peaks)))
{
  MassTrace tmp_mt(peak_vec);

  std::vector<PeakType>::const_iterator v_it = peak_vec.begin();

  for (MassTrace::const_iterator m_it = tmp_mt.begin(); m_it != tmp_mt.end(); ++m_it)
  {
      TEST_EQUAL(*v_it, *m_it);
      ++v_it;
  }

  TEST_REAL_SIMILAR(tmp_mt.getAverageMS1CycleTime(), (tmp_peak6.getRT() - tmp_peak0.getRT()) / (7-1));
}
END_SECTION

/////



MassTrace test_mt(peak_lst);
test_mt.updateWeightedMeanRT();
test_mt.updateWeightedMeanMZ();


/////////////////////////////////////////////////////////////
// operator tests
/////////////////////////////////////////////////////////////

START_SECTION((PeakType& operator[](const Size &mt_idx)))
{
  TEST_REAL_SIMILAR(test_mt[1].getRT(), 153.23);
  TEST_REAL_SIMILAR(test_mt[1].getMZ(), 230.10235);
  TEST_REAL_SIMILAR(test_mt[1].getIntensity(), (OpenMS::Peak2D::IntensityType)542293.0);

  TEST_REAL_SIMILAR(test_mt[4].getRT(), 156.233);
  TEST_REAL_SIMILAR(test_mt[4].getMZ(), 230.10116);
  TEST_REAL_SIMILAR(test_mt[4].getIntensity(), (OpenMS::Peak2D::IntensityType)17342933.0);
}
END_SECTION

/////

START_SECTION((const PeakType& operator[](const Size &mt_idx) const ))
{
  const MassTrace test_mt_const(test_mt);

  double rt1 = test_mt_const[1].getRT();
  double mz1 = test_mt_const[1].getMZ();
  double int1 = test_mt_const[1].getIntensity();
  double rt2 = test_mt_const[4].getRT();
  double mz2 = test_mt_const[4].getMZ();
  double int2 = test_mt_const[4].getIntensity();

  TEST_REAL_SIMILAR(rt1, 153.23);
  TEST_REAL_SIMILAR(mz1, 230.10235);
  TEST_REAL_SIMILAR(int1, 542293.0);

  TEST_REAL_SIMILAR(rt2, 156.233);
  TEST_REAL_SIMILAR(mz2, 230.10116);
  TEST_REAL_SIMILAR(int2, 17342933.0);
}
END_SECTION


/////////////////////////////////////////////////////////////
// iterator tests
/////////////////////////////////////////////////////////////

std::vector<PeakType>::iterator start_it = peak_vec.begin();
std::vector<PeakType>::iterator end_it = peak_vec.end();
std::vector<PeakType>::reverse_iterator rstart_it = peak_vec.rbegin();
std::vector<PeakType>::reverse_iterator rend_it = peak_vec.rend();

std::vector<PeakType>::const_iterator const_start_it = peak_vec.begin();
std::vector<PeakType>::const_iterator const_end_it = peak_vec.end();
std::vector<PeakType>::const_reverse_iterator const_rstart_it = peak_vec.rbegin();
std::vector<PeakType>::const_reverse_iterator const_rend_it = peak_vec.rend();


START_SECTION((iterator begin()))
{
    MassTrace::iterator mt_it = test_mt.begin();
    TEST_EQUAL(*start_it, *mt_it);
}
END_SECTION

/////

START_SECTION((iterator end()))
{
  MassTrace::iterator mt_it = test_mt.begin();

  for ( ; mt_it != test_mt.end() - 1; ++mt_it)
  {
  }
  TEST_EQUAL(*(end_it - 1), *mt_it);
}
END_SECTION

/////

START_SECTION((const_iterator begin() const))
{
  MassTrace::const_iterator mt_it = test_mt.begin();
  TEST_EQUAL(*const_start_it, *mt_it);
}
END_SECTION

/////

START_SECTION((const_iterator end() const))
{
  MassTrace::const_iterator mt_it = test_mt.begin();

  for ( ; mt_it != test_mt.end() - 1; ++mt_it)
  {
  }
  TEST_EQUAL(*(const_end_it - 1), *mt_it);
}
END_SECTION

/////

START_SECTION((reverse_iterator rbegin()))
{
  MassTrace::reverse_iterator mt_it = test_mt.rbegin();
  TEST_EQUAL(*rstart_it, *mt_it);
}
END_SECTION

/////

START_SECTION((reverse_iterator rend()))
{
  MassTrace::reverse_iterator mt_it = test_mt.rbegin();

  for ( ; mt_it != test_mt.rend() - 1; ++mt_it)
  {
  }
  TEST_EQUAL(*(rend_it - 1), *mt_it);
}
END_SECTION

/////

START_SECTION((const_reverse_iterator rbegin() const))
{
  MassTrace::const_reverse_iterator mt_it = test_mt.rbegin();
  TEST_EQUAL(*const_rstart_it, *mt_it);
}
END_SECTION

/////

START_SECTION((const_reverse_iterator rend() const))
{
  MassTrace::const_reverse_iterator mt_it = test_mt.rbegin();

  for ( ; mt_it != test_mt.rend() - 1; ++mt_it)
  {
  }
  TEST_EQUAL(*(const_rend_it - 1), *mt_it);
}
END_SECTION

/////////////////////////////////////////////////////////////
// accessor method tests
/////////////////////////////////////////////////////////////


START_SECTION((Size getSize() const))
{
  Size test_mt_size = test_mt.getSize();
  TEST_EQUAL(test_mt_size, 7);
}
END_SECTION

/////

START_SECTION((String getLabel() const ))
{
  const String test_mt_label = test_mt.getLabel();
  TEST_EQUAL(test_mt_label, "");
}
END_SECTION

/////

START_SECTION((void setLabel(const String& label)))
{
  test_mt.setLabel("TEST_TRACE");
  String test_mt_label = test_mt.getLabel();

  TEST_EQUAL(test_mt_label, "TEST_TRACE");
}
END_SECTION

/////


START_SECTION((double getCentroidMZ() const ))
{
  MassTrace test_mt_const(test_mt);
  double test_mt_cent_mz = test_mt_const.getCentroidMZ();
  TEST_REAL_SIMILAR(test_mt_cent_mz, 230.10188);
}
END_SECTION

/////
START_SECTION((double getCentroidRT() const ))
{
  MassTrace test_mt_const(test_mt);
  double test_mt_cent_rt = test_mt_const.getCentroidRT();
  TEST_REAL_SIMILAR(test_mt_cent_rt, 155.214671250425);
}
END_SECTION

/////

START_SECTION((double getAverageMS1CycleTime() const))
{
  MassTrace tmp_mt(peak_lst);
  TEST_REAL_SIMILAR(tmp_mt.getAverageMS1CycleTime(), (tmp_peak6.getRT() - tmp_peak0.getRT()) / (7-1));
}
END_SECTION

/////


START_SECTION((void updateWeightedMZsd()))
{
  MassTrace empty_trace;
  TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateWeightedMZsd());
  std::vector<PeakType> peaks;
  PeakType p1, p2;
  p1.setMZ(123.123);
  p1.setIntensity(0.0);
  p2.setMZ(123.321);
  p2.setIntensity(0.0);
  peaks.push_back(p1);
  peaks.push_back(p2);
  MassTrace zero_int_mt(peaks);
  TEST_EXCEPTION(Exception::InvalidValue, zero_int_mt.updateWeightedMZsd());

  test_mt.updateWeightedMZsd();
  double test_mt_sd = test_mt.getCentroidSD();

  TEST_REAL_SIMILAR(test_mt_sd, 0.0004594);
}
END_SECTION

/////
                          
START_SECTION((double getCentroidSD() const ))
{
  MassTrace test_mt_const(test_mt);
  double test_mt_sd = test_mt_const.getCentroidSD();
  TEST_REAL_SIMILAR(test_mt_sd, 0.0004594);
}
END_SECTION


/////

START_SECTION((double getTraceLength() const ))
{
  const MassTrace test_mt_const(test_mt);

  double mt_length = test_mt_const.getTraceLength();

  TEST_REAL_SIMILAR(mt_length, 6.018)
}
END_SECTION

/////

std::vector<double> smoothed_ints;
smoothed_ints.push_back(500.0);
smoothed_ints.push_back(540000.0);
smoothed_ints.push_back(18000000.0);
smoothed_ints.push_back(33000000.0);
smoothed_ints.push_back(17500000.0);
smoothed_ints.push_back(540000.0);
smoothed_ints.push_back(549223.0);
smoothed_ints.push_back(300.0);


START_SECTION((void setSmoothedIntensities(const std::vector<double>& db_vec)))
{
  TEST_EXCEPTION(Exception::InvalidValue, test_mt.setSmoothedIntensities(smoothed_ints));
  smoothed_ints.pop_back();

  test_mt.setSmoothedIntensities(smoothed_ints);
  TEST_EQUAL(test_mt.getSmoothedIntensities().size(), smoothed_ints.size())
}
END_SECTION

/////

START_SECTION((std::vector<double> getSmoothedIntensities()))
{
  std::vector<double> smoothed_vec = test_mt.getSmoothedIntensities();

  TEST_EQUAL(smoothed_vec.empty(), false);
  TEST_EQUAL(smoothed_vec.size(), smoothed_ints.size());
}
END_SECTION

/////

test_mt.setSmoothedIntensities(smoothed_ints);

START_SECTION((double getIntensity(bool smoothed) const))
{
  test_mt.estimateFWHM(true);

  double smoothed_area = test_mt.getIntensity(true);
  TEST_REAL_SIMILAR(smoothed_area, 69505990.0000001);

  double raw_area = test_mt.getIntensity(false);
  TEST_REAL_SIMILAR(raw_area, 69863097.2125001);
}
END_SECTION

/////

START_SECTION((double getMaxIntensity(bool smoothed) const))
{
  double smoothed_maxint = test_mt.getMaxIntensity(true);
  TEST_REAL_SIMILAR(smoothed_maxint, 33000000.0);

  double raw_maxint= test_mt.getMaxIntensity(false);
  TEST_REAL_SIMILAR(raw_maxint, 33329536.0);
}
END_SECTION

/////

START_SECTION((double getMaxIntensity(bool) const ))
{
  const MassTrace test_mt_const(test_mt);
  double smoothed_maxint = test_mt_const.getMaxIntensity(true);
  TEST_REAL_SIMILAR(smoothed_maxint, 33000000.0);

  double raw_maxint= test_mt_const.getMaxIntensity(false);
  TEST_REAL_SIMILAR(raw_maxint, 33329536.0);
}
END_SECTION

/////

START_SECTION((const std::vector<double>& getSmoothedIntensities() const))
{
  std::vector<double> smoothed_vec = test_mt.getSmoothedIntensities();

  TEST_EQUAL(smoothed_vec.empty(), false);
  TEST_EQUAL(smoothed_vec.size(), smoothed_ints.size());
}
END_SECTION

/////

MassTrace test_mt2(peak_vec), test_mt3;
test_mt2.updateWeightedMeanRT();
test_mt2.updateWeightedMeanMZ();


START_SECTION((double getFWHM() const))
{
  double test_mt_fwhm = test_mt.getFWHM();
  TEST_REAL_SIMILAR(test_mt_fwhm, 3.05481743986255);
}
END_SECTION
                                         

/////

START_SECTION((double computeSmoothedPeakArea() const))
{
  double peak_area = test_mt.computeSmoothedPeakArea();
  TEST_REAL_SIMILAR(peak_area, 70303689.0475001)
}
END_SECTION

/////

// create Masstraces which can be used to test missing RT scans
std::vector<PeakType> peak_vec_red;
double rt[] = {152.0,  153.0,       154.0,         155.0,        156.0,        157.0,        158.0,  159.0, 160.0};
double it[] = {  542, 542293, 1.82824e+007, 3.33295e+007, 3.33295e+007, 3.33295e+007, 1.73429e+007, 333291,   339};
for (Size i=0; i<9; ++i)
{
  peak_vec_red.push_back(fillPeak(rt[i], 100.0, it[i]));
}
std::vector<PeakType> peak_vec_red2 = peak_vec_red;
peak_vec_red2.erase(peak_vec_red2.begin() + 4);
MassTrace tmp_mt(peak_vec_red);
MassTrace tmp_mt2(peak_vec_red2);

START_SECTION((double computePeakArea() const))
{
  double peak_area = test_mt.computePeakArea();
  TEST_REAL_SIMILAR(peak_area, 70303710.2575001)

  // ensure that missing peaks (=scans) do not impact quantification
  TEST_REAL_SIMILAR(tmp_mt.computePeakArea(), tmp_mt2.computePeakArea())
}
END_SECTION

/////

START_SECTION((double computeFwhmAreaSmooth() const))
{
  double peak_area = test_mt.computeFwhmAreaSmooth();
  TEST_REAL_SIMILAR(peak_area, 69505990.0000001)
}
END_SECTION

/////

START_SECTION((double computeFwhmArea() const))
{
  double peak_area = test_mt.computeFwhmArea();
  TEST_REAL_SIMILAR(peak_area, 69863097.2125001)

  // ensure that missing peaks (=scans) do not impact quantification
  tmp_mt.estimateFWHM();
  tmp_mt2.estimateFWHM();
  TEST_REAL_SIMILAR(tmp_mt.computeFwhmArea(), tmp_mt2.computeFwhmArea())
}
END_SECTION


/////

START_SECTION((Size findMaxByIntPeak(bool use_smoothed_ints = false) const ))
{
  TEST_EXCEPTION(Exception::InvalidValue, test_mt2.findMaxByIntPeak(true));
  TEST_EXCEPTION(Exception::InvalidValue, test_mt3.findMaxByIntPeak(false));
  TEST_EXCEPTION(Exception::InvalidValue, test_mt3.findMaxByIntPeak(true));

  Size max_peak_idx1 = test_mt.findMaxByIntPeak(true);
  Size max_peak_idx2 = test_mt.findMaxByIntPeak(false);

  TEST_EQUAL(max_peak_idx1, 3);
  TEST_EQUAL(max_peak_idx2, 3);
}
END_SECTION

/////

START_SECTION((double estimateFWHM(bool use_smoothed_ints = false)))
{
  TEST_EXCEPTION(Exception::InvalidValue, test_mt2.estimateFWHM(true));
  TEST_EXCEPTION(Exception::InvalidValue, test_mt3.estimateFWHM(false));

  double test_fwhm1 = test_mt.estimateFWHM(false);
  double test_fwhm2 = test_mt.estimateFWHM(true);

  TEST_REAL_SIMILAR(test_fwhm1, 3.07921244222942);
  TEST_REAL_SIMILAR(test_fwhm2, 3.05481743986255);
}
END_SECTION

/////
START_SECTION(static MT_QUANTMETHOD getQuantMethod(const String& val))
  
  TEST_EQUAL(MassTrace::getQuantMethod("area"), MassTrace::MT_QUANT_AREA)
  TEST_EQUAL(MassTrace::getQuantMethod("median"), MassTrace::MT_QUANT_MEDIAN)
  TEST_EQUAL(MassTrace::getQuantMethod("max_height"), MassTrace::MT_QUANT_HEIGHT)
  TEST_EQUAL(MassTrace::getQuantMethod("somethingwrong"), MassTrace::SIZE_OF_MT_QUANTMETHOD)

END_SECTION

START_SECTION((void setQuantMethod(MT_QUANTMETHOD method)))
{
  MassTrace mt_empty;
  TEST_EQUAL(mt_empty.getQuantMethod(), MassTrace::MT_QUANT_AREA);

  MassTrace raw_mt(peak_vec);
  // area (default)
  raw_mt.estimateFWHM(false);
  TEST_REAL_SIMILAR(raw_mt.getIntensity(false), 69863097.2125001);

  raw_mt.setQuantMethod(MassTrace::MT_QUANT_MEDIAN);
  // should return the median of the intensities
  TEST_REAL_SIMILAR(raw_mt.getIntensity(false), 542293.0);
  TEST_EQUAL(raw_mt.getQuantMethod(), MassTrace::MT_QUANT_MEDIAN);

  // test max_height
  raw_mt.setQuantMethod(MassTrace::MT_QUANT_HEIGHT);
  TEST_REAL_SIMILAR(raw_mt.getIntensity(false), 33329535.0);
  TEST_EQUAL(raw_mt.getQuantMethod(), MassTrace::MT_QUANT_HEIGHT);

  TEST_EXCEPTION(Exception::InvalidValue, raw_mt.setQuantMethod(MassTrace::SIZE_OF_MT_QUANTMETHOD))

}
END_SECTION

START_SECTION((MT_QUANTMETHOD getQuantMethod() const))
{
  NOT_TESTABLE // tested above
}
END_SECTION
/////

START_SECTION((std::pair<Size, Size> getFWHMborders() const ))
{
  MassTrace raw_mt(peak_vec);
  std::pair<Size, Size> interval = raw_mt.getFWHMborders();

  TEST_EQUAL(interval.first, 0);
  TEST_EQUAL(interval.second, 0);

  interval = test_mt.getFWHMborders();

  TEST_EQUAL(interval.first, 1);
  TEST_EQUAL(interval.second, 5);
}
END_SECTION

/////

std::vector<PeakType> double_peak(peak_vec);
double_peak.insert(double_peak.end(), peak_vec.begin(), peak_vec.end());

std::vector<double> double_smooth_ints(smoothed_ints);
double_smooth_ints.insert(double_smooth_ints.end(), smoothed_ints.begin(), smoothed_ints.end());

MassTrace double_mt(double_peak);
double_mt.setSmoothedIntensities(double_smooth_ints);


START_SECTION((MassTrace(const MassTrace &)))
{
  MassTrace copy_mt(test_mt);

  MassTrace::const_iterator c_it = copy_mt.begin();

  for (MassTrace::const_iterator t_it = test_mt.begin(); t_it != test_mt.end(); ++t_it)
  {
    TEST_EQUAL(*c_it, *t_it);
    ++c_it;
  }

  TEST_REAL_SIMILAR(copy_mt.getCentroidMZ(), test_mt.getCentroidMZ());
  TEST_REAL_SIMILAR(copy_mt.getCentroidRT(), test_mt.getCentroidRT());

  TEST_EQUAL(copy_mt.getLabel(), test_mt.getLabel());

  std::vector<double> sm1(copy_mt.getSmoothedIntensities()), sm2(test_mt.getSmoothedIntensities());

  std::vector<double>::const_iterator sm1_it = sm1.begin();

  for (std::vector<double>::const_iterator sm2_it = sm2.begin(); sm2_it != sm2.end(); ++sm2_it)
  {
    TEST_EQUAL(*sm1_it, *sm2_it);
    ++sm1_it;
  }
}
END_SECTION

/////

START_SECTION((MassTrace& operator=(const MassTrace &)))
{
  MassTrace copy_mt = test_mt;

  MassTrace::const_iterator c_it = copy_mt.begin();

  for (MassTrace::const_iterator t_it = test_mt.begin(); t_it != test_mt.end(); ++t_it)
  {
      TEST_EQUAL(*c_it, *t_it);
      ++c_it;
  }

  TEST_REAL_SIMILAR(copy_mt.getCentroidMZ(), test_mt.getCentroidMZ());
  TEST_REAL_SIMILAR(copy_mt.getCentroidRT(), test_mt.getCentroidRT());
  TEST_EQUAL(copy_mt.getLabel(), test_mt.getLabel());

  std::vector<double> sm1(copy_mt.getSmoothedIntensities()), sm2(test_mt.getSmoothedIntensities());
  std::vector<double>::const_iterator sm1_it = sm1.begin();
  for (std::vector<double>::const_iterator sm2_it = sm2.begin(); sm2_it != sm2.end(); ++sm2_it)
  {
    TEST_EQUAL(*sm1_it, *sm2_it);
    ++sm1_it;
  }
}
END_SECTION

/////

START_SECTION((ConvexHull2D getConvexhull() const ))
{
  ConvexHull2D tmp_hull = test_mt.getConvexhull();
  DPosition<2> tmp_p1(154.21, 230.10181);
  DPosition<2> tmp_p2(155.22, 230.10181);
  DPosition<2> tmp_p3(154.21, 229.10181);

  TEST_EQUAL(tmp_hull.encloses(tmp_p1), true);
  TEST_EQUAL(tmp_hull.encloses(tmp_p2), false);
  TEST_EQUAL(tmp_hull.encloses(tmp_p3), false);
}
END_SECTION

/////

MassTrace empty_trace;

START_SECTION((void updateWeightedMeanRT()))
{
  TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateWeightedMeanRT());

  test_mt.updateWeightedMeanRT();
  TEST_REAL_SIMILAR(test_mt.getCentroidRT(), 155.214671250425);
}
END_SECTION

/////

START_SECTION((void updateMedianRT()))
{
  TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateMedianRT());

  test_mt.updateMedianRT();
  TEST_REAL_SIMILAR(test_mt.getCentroidRT(), 155.24);
}
END_SECTION

/////

START_SECTION((void updateMedianMZ()))
{
  TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateMedianMZ());

  test_mt.updateMedianMZ();
  TEST_REAL_SIMILAR(test_mt.getCentroidMZ(), 230.10198);
}
END_SECTION

/////

START_SECTION((void updateMeanMZ()))
{
  TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateMeanMZ());

  test_mt.updateMeanMZ();
  TEST_REAL_SIMILAR(test_mt.getCentroidMZ(), 230.101918);
}
END_SECTION

/////

START_SECTION((void updateWeightedMeanMZ()))
{
  TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateWeightedMeanMZ());

  test_mt.updateWeightedMeanMZ();
  TEST_REAL_SIMILAR(test_mt.getCentroidMZ(), 230.101883054967);
}
END_SECTION

/////

START_SECTION((void updateSmoothedMaxRT()))
{
  MassTrace raw_mt(peak_vec);
  TEST_EXCEPTION(Exception::InvalidValue, raw_mt.updateSmoothedMaxRT());

  test_mt.updateSmoothedMaxRT();
  double smooth_max_rt = test_mt.getCentroidRT();
  TEST_REAL_SIMILAR(smooth_max_rt, 155.24);
}
END_SECTION

/////

START_SECTION((void updateSmoothedWeightedMeanRT()))
{
    MassTrace raw_mt(peak_vec);
    TEST_EXCEPTION(Exception::InvalidValue, raw_mt.updateSmoothedWeightedMeanRT());

    test_mt.updateSmoothedWeightedMeanRT();
    double smooth_max_rt = test_mt.getCentroidRT();
    TEST_REAL_SIMILAR(smooth_max_rt, 155.2468039);
}
END_SECTION

/////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
