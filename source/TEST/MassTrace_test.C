// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MassTrace.h>

///////////////////////////

START_TEST(MassTrace, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MassTrace* d_ptr = 0;
MassTrace* nullPointer = 0;
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

PeakType tmp_peak1, tmp_peak2, tmp_peak3, tmp_peak4, tmp_peak5;

tmp_peak1.setRT(153.23);
tmp_peak1.setMZ(230.10235);
tmp_peak1.setIntensity(54392293.0);
peak_vec.push_back(tmp_peak1);
peak_lst.push_back(tmp_peak1);

tmp_peak2.setRT(154.21);
tmp_peak2.setMZ(230.10181);
tmp_peak2.setIntensity(1828482393.0);
peak_vec.push_back(tmp_peak2);
peak_lst.push_back(tmp_peak2);

tmp_peak3.setRT(155.24);
tmp_peak3.setMZ(230.10229);
tmp_peak3.setIntensity(3339229535.0);
peak_vec.push_back(tmp_peak3);
peak_lst.push_back(tmp_peak3);

tmp_peak4.setRT(156.233);
tmp_peak4.setMZ(230.10116);
tmp_peak4.setIntensity(1734222933.0);
peak_vec.push_back(tmp_peak4);
peak_lst.push_back(tmp_peak4);

tmp_peak5.setRT(157.24);
tmp_peak5.setMZ(230.10198);
tmp_peak5.setIntensity(33392291.0);
peak_vec.push_back(tmp_peak5);
peak_lst.push_back(tmp_peak5);


/////////////////////////////////////////////////////////////
// detailed constructors test
/////////////////////////////////////////////////////////////

START_SECTION((MassTrace(const std::list< PeakType > &)))
{
    MassTrace tmp_mt(peak_lst);

    std::list<PeakType>::const_iterator l_it = peak_lst.begin();

    for (MassTrace::const_iterator m_it = tmp_mt.begin(); m_it != tmp_mt.end(); ++m_it)
    {
        TEST_EQUAL(*l_it, *m_it);
        ++l_it;
    }
}
END_SECTION

/////

START_SECTION((MassTrace(const std::vector< PeakType > &)))
{
    MassTrace tmp_mt(peak_vec);

    std::vector<PeakType>::const_iterator v_it = peak_vec.begin();

    for (MassTrace::const_iterator m_it = tmp_mt.begin(); m_it != tmp_mt.end(); ++m_it)
    {
        TEST_EQUAL(*v_it, *m_it);
        ++v_it;
    }
}
END_SECTION

/////



MassTrace test_mt(peak_lst);
test_mt.updateWeightedMeanRT();
test_mt.updateWeightedMeanMZ();



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

    TEST_EQUAL(test_mt_size, 5);
}
END_SECTION

/////

START_SECTION((String getLabel()))
{
    String test_mt_label = test_mt.getLabel();

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

// std::cerr << setprecision(10) << test_mt_cent_mz << std::endl;

START_SECTION((DoubleReal getCentroidMZ()))
{
    DoubleReal test_mt_cent_mz = test_mt.getCentroidMZ();

    TEST_REAL_SIMILAR(test_mt_cent_mz, 230.1022624);
}
END_SECTION

/////

START_SECTION((DoubleReal getCentroidRT()))
{
    DoubleReal test_mt_cent_rt = test_mt.getCentroidRT();

    std::cerr << setprecision(10) << test_mt_cent_rt << std::endl;

    TEST_REAL_SIMILAR(test_mt_cent_rt, 155.2108433);
}
END_SECTION

/////

std::vector<DoubleReal> smoothed_ints;
smoothed_ints.push_back(54000000.0);
smoothed_ints.push_back(1800000000.0);
smoothed_ints.push_back(3300000000.0);
smoothed_ints.push_back(1750000000.0);
smoothed_ints.push_back(54000000.0);
smoothed_ints.push_back(54922953.0);

START_SECTION((void setSmoothedIntensities(const std::vector<DoubleReal>& db_vec)))
{
    TEST_EXCEPTION(Exception::InvalidValue, test_mt.setSmoothedIntensities(smoothed_ints));
    smoothed_ints.pop_back();

    test_mt.setSmoothedIntensities(smoothed_ints);
}
END_SECTION

/////

START_SECTION((std::vector<DoubleReal> getSmoothedIntensities()))
{
    std::vector<DoubleReal> smoothed_vec = test_mt.getSmoothedIntensities();

    TEST_EQUAL(smoothed_vec.empty(), false);
    TEST_EQUAL(smoothed_vec.size(), smoothed_ints.size());
}
END_SECTION

/////

START_SECTION((std::vector<DoubleReal> getSmoothedIntensities() const))
{
    std::vector<DoubleReal> smoothed_vec = test_mt.getSmoothedIntensities();

    TEST_EQUAL(smoothed_vec.empty(), false);
    TEST_EQUAL(smoothed_vec.size(), smoothed_ints.size());
}
END_SECTION

/////

MassTrace test_mt2(peak_vec), test_mt3;
test_mt2.updateWeightedMeanRT();
test_mt2.updateWeightedMeanMZ();


START_SECTION((void setFWHMScansNum(Size r_fwhm)))
{
    test_mt.setFWHMScansNum(2);
    TEST_EQUAL(test_mt.getFWHMScansNum(), 2);
}
END_SECTION

/////

START_SECTION((Size getFWHMScansNum()))
{
    NOT_TESTABLE; // see setFWHMScansNum
}
END_SECTION

/////

START_SECTION((DoubleReal computePeakArea()))
{
    DoubleReal peak_area = test_mt.computePeakArea();

    TEST_REAL_SIMILAR(peak_area, 6989719445.0)
}
END_SECTION

/////

START_SECTION((Size findMaxByIntPeak(bool)))
{
    TEST_EXCEPTION(Exception::InvalidValue, test_mt2.findMaxByIntPeak(true));
    TEST_EXCEPTION(Exception::InvalidValue, test_mt3.findMaxByIntPeak(false));
    TEST_EXCEPTION(Exception::InvalidValue, test_mt3.findMaxByIntPeak(true));

    Size max_peak_idx1 = test_mt.findMaxByIntPeak(true);
    Size max_peak_idx2 = test_mt.findMaxByIntPeak(false);

    TEST_EQUAL(max_peak_idx1, 2);
    TEST_EQUAL(max_peak_idx2, 2);
}
END_SECTION

/////

START_SECTION((DoubleReal estimateFWHM(bool)))
{
    TEST_EXCEPTION(Exception::InvalidValue, test_mt2.estimateFWHM(true));
    TEST_EXCEPTION(Exception::InvalidValue, test_mt3.estimateFWHM(false));

    DoubleReal test_fwhm1 = test_mt.estimateFWHM(false);
    DoubleReal test_fwhm2 = test_mt.estimateFWHM(true);

    TEST_REAL_SIMILAR(test_fwhm1, 4.01);
    TEST_REAL_SIMILAR(test_fwhm2, 4.01);

    TEST_EQUAL(test_mt.getFWHMScansNum(), 5);
}
END_SECTION

/////

std::vector<PeakType> double_peak(peak_vec);
double_peak.insert(double_peak.end(), peak_vec.begin(), peak_vec.end());

std::vector<DoubleReal> double_smooth_ints(smoothed_ints);
double_smooth_ints.insert(double_smooth_ints.end(), smoothed_ints.begin(), smoothed_ints.end());

MassTrace double_mt(double_peak);
double_mt.setSmoothedIntensities(double_smooth_ints);


START_SECTION((void findLocalExtrema(const Size &, std::vector< Size > &, std::vector< Size > &)))
{
    std::vector<Size> maxes, mins;

    double_mt.findLocalExtrema(2, maxes, mins);

    TEST_EQUAL(maxes.size(), 2);
    TEST_EQUAL(mins.size(), 1);

    maxes.clear();
    mins.clear();

    double_mt.findLocalExtrema(3, maxes, mins);

    TEST_EQUAL(maxes.size(), 1);
    TEST_EQUAL(mins.size(), 0);

    // std::cerr << maxes.size() << " " << mins.size() << std::endl;
}
END_SECTION


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

    std::vector<DoubleReal> sm1(copy_mt.getSmoothedIntensities()), sm2(test_mt.getSmoothedIntensities());

    std::vector<DoubleReal>::const_iterator sm1_it = sm1.begin();

    for (std::vector<DoubleReal>::const_iterator sm2_it = sm2.begin(); sm2_it != sm2.end(); ++sm2_it)
    {
        TEST_EQUAL(*sm1_it, *sm2_it);
        ++sm1_it;
    }


    TEST_EQUAL(copy_mt.getFWHMScansNum(), test_mt.getFWHMScansNum());
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

    std::vector<DoubleReal> sm1(copy_mt.getSmoothedIntensities()), sm2(test_mt.getSmoothedIntensities());

    std::vector<DoubleReal>::const_iterator sm1_it = sm1.begin();

    for (std::vector<DoubleReal>::const_iterator sm2_it = sm2.begin(); sm2_it != sm2.end(); ++sm2_it)
    {
        TEST_EQUAL(*sm1_it, *sm2_it);
        ++sm1_it;
    }


    TEST_EQUAL(copy_mt.getFWHMScansNum(), test_mt.getFWHMScansNum());
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

    TEST_REAL_SIMILAR(test_mt.getCentroidRT(),155.210843262781);
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
