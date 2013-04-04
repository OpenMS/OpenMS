// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

PeakType tmp_peak0, tmp_peak1, tmp_peak2, tmp_peak3, tmp_peak4, tmp_peak5, tmp_peak6;

tmp_peak0.setRT(152.22);
tmp_peak0.setMZ(230.10223);
tmp_peak0.setIntensity(542.0);
peak_vec.push_back(tmp_peak0);
peak_lst.push_back(tmp_peak0);

tmp_peak1.setRT(153.23);
tmp_peak1.setMZ(230.10235);
tmp_peak1.setIntensity(542293.0);
peak_vec.push_back(tmp_peak1);
peak_lst.push_back(tmp_peak1);

tmp_peak2.setRT(154.21);
tmp_peak2.setMZ(230.10181);
tmp_peak2.setIntensity(18282393.0);
peak_vec.push_back(tmp_peak2);
peak_lst.push_back(tmp_peak2);

tmp_peak3.setRT(155.24);
tmp_peak3.setMZ(230.10229);
tmp_peak3.setIntensity(33329535.0);
peak_vec.push_back(tmp_peak3);
peak_lst.push_back(tmp_peak3);

tmp_peak4.setRT(156.233);
tmp_peak4.setMZ(230.10116);
tmp_peak4.setIntensity(17342933.0);
peak_vec.push_back(tmp_peak4);
peak_lst.push_back(tmp_peak4);

tmp_peak5.setRT(157.24);
tmp_peak5.setMZ(230.10198);
tmp_peak5.setIntensity(333291.0);
peak_vec.push_back(tmp_peak5);
peak_lst.push_back(tmp_peak5);

tmp_peak6.setRT(158.238);
tmp_peak6.setMZ(230.10254);
tmp_peak6.setIntensity(339.0);
peak_vec.push_back(tmp_peak5);
peak_lst.push_back(tmp_peak5);


/////////////////////////////////////////////////////////////
// detailed constructors test
/////////////////////////////////////////////////////////////

START_SECTION((MassTrace(const std::list< PeakType > &, const DoubleReal &scan_time=1.0)))
{
    MassTrace tmp_mt(peak_lst);

    std::list<PeakType>::const_iterator l_it = peak_lst.begin();

    for (MassTrace::const_iterator m_it = tmp_mt.begin(); m_it != tmp_mt.end(); ++m_it)
    {
        TEST_EQUAL(*l_it, *m_it);
        ++l_it;
    }

    TEST_REAL_SIMILAR(tmp_mt.getScanTime(), 1.0);

    MassTrace tmp_mt2(peak_lst, 0.25);
    TEST_REAL_SIMILAR(tmp_mt2.getScanTime(), 0.25);


}
END_SECTION

/////

START_SECTION((MassTrace(const std::vector< PeakType > &, const DoubleReal &scan_time=1.0)))
{
    MassTrace tmp_mt(peak_vec);

    std::vector<PeakType>::const_iterator v_it = peak_vec.begin();

    for (MassTrace::const_iterator m_it = tmp_mt.begin(); m_it != tmp_mt.end(); ++m_it)
    {
        TEST_EQUAL(*v_it, *m_it);
        ++v_it;
    }

    TEST_REAL_SIMILAR(tmp_mt.getScanTime(), 1.0);

    MassTrace tmp_mt2(peak_lst, 0.25);
    TEST_REAL_SIMILAR(tmp_mt2.getScanTime(), 0.25);
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
    TEST_REAL_SIMILAR(test_mt[1].getIntensity(), 542293.0);

    TEST_REAL_SIMILAR(test_mt[4].getRT(), 156.233);
    TEST_REAL_SIMILAR(test_mt[4].getMZ(), 230.10116);
    TEST_REAL_SIMILAR(test_mt[4].getIntensity(), 17342933.0);
}
END_SECTION

/////

START_SECTION((const PeakType& operator[](const Size &mt_idx) const ))
{
    const MassTrace test_mt_const(test_mt);

    DoubleReal rt1 = test_mt_const[1].getRT();
    DoubleReal mz1 = test_mt_const[1].getMZ();
    DoubleReal int1 = test_mt_const[1].getIntensity();
    DoubleReal rt2 = test_mt_const[4].getRT();
    DoubleReal mz2 = test_mt_const[4].getMZ();
    DoubleReal int2 = test_mt_const[4].getIntensity();

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



START_SECTION((DoubleReal getCentroidMZ()))
{
    DoubleReal test_mt_cent_mz = test_mt.getCentroidMZ();

    TEST_REAL_SIMILAR(test_mt_cent_mz, 230.10188);
}
END_SECTION

/////

START_SECTION((DoubleReal getCentroidMZ() const ))
{
    MassTrace test_mt_const(test_mt);

    DoubleReal test_mt_cent_mz = test_mt_const.getCentroidMZ();

    TEST_REAL_SIMILAR(test_mt_cent_mz, 230.10188);
}
END_SECTION

/////
START_SECTION((DoubleReal getCentroidRT()))
{
    DoubleReal test_mt_cent_rt = test_mt.getCentroidRT();

    TEST_REAL_SIMILAR(test_mt_cent_rt, 155.2205);
}
END_SECTION

/////
START_SECTION((DoubleReal getCentroidRT() const ))
{
    MassTrace test_mt_const(test_mt);

    DoubleReal test_mt_cent_rt = test_mt_const.getCentroidRT();

    TEST_REAL_SIMILAR(test_mt_cent_rt, 155.2205);
}
END_SECTION

/////

START_SECTION((DoubleReal getScanTime()))
{
    MassTrace tmp_mt(peak_lst, 0.25);

    DoubleReal test_scantime = tmp_mt.getScanTime();

    TEST_REAL_SIMILAR(test_scantime, 0.25);
}
END_SECTION

/////



PeakType p1, p2;
p1.setMZ(123.123);
p1.setIntensity(0.0);
p2.setMZ(123.321);
p2.setIntensity(0.0);

std::vector<PeakType> peaks;
peaks.push_back(p1);
peaks.push_back(p2);

MassTrace zero_int_mt(peaks);

START_SECTION((void updateWeightedMZsd()))
{
    MassTrace empty_trace;
    TEST_EXCEPTION(Exception::InvalidValue, empty_trace.updateWeightedMZsd());

    test_mt.updateWeightedMZsd();
    DoubleReal test_mt_sd = test_mt.getCentroidSD();

    TEST_REAL_SIMILAR(test_mt_sd, 0.0004594);

    TEST_EXCEPTION(Exception::InvalidValue, zero_int_mt.updateWeightedMZsd());
}
END_SECTION

/////

START_SECTION((DoubleReal getCentroidSD()))
{
    DoubleReal test_mt_sd = test_mt.getCentroidSD();

    TEST_REAL_SIMILAR(test_mt_sd, 0.0004594);
}
END_SECTION

/////

START_SECTION((DoubleReal getCentroidSD() const ))
{
    MassTrace test_mt_const(test_mt);

    DoubleReal test_mt_sd = test_mt_const.getCentroidSD();

    TEST_REAL_SIMILAR(test_mt_sd, 0.0004594);
}
END_SECTION

/////
START_SECTION((void setCentroidSD(const DoubleReal &tmp_sd)))
{
    test_mt.setCentroidSD(0.00048);
    DoubleReal test_mt_sd = test_mt.getCentroidSD();

    TEST_REAL_SIMILAR(test_mt_sd, 0.00048);
}
END_SECTION

/////

START_SECTION((DoubleReal getTraceLength()))
{
    DoubleReal mt_length = test_mt.getTraceLength();

    TEST_REAL_SIMILAR(mt_length, 5.02)
}
END_SECTION

/////

START_SECTION((DoubleReal getTraceLength() const ))
{
    MassTrace test_mt_const(test_mt);

    DoubleReal mt_length = test_mt_const.getTraceLength();

    TEST_REAL_SIMILAR(mt_length, 5.02)
}
END_SECTION

/////

std::vector<DoubleReal> smoothed_ints;
smoothed_ints.push_back(500.0);
smoothed_ints.push_back(540000.0);
smoothed_ints.push_back(18000000.0);
smoothed_ints.push_back(33000000.0);
smoothed_ints.push_back(17500000.0);
smoothed_ints.push_back(540000.0);
smoothed_ints.push_back(549223.0);
smoothed_ints.push_back(300.0);


START_SECTION((void setSmoothedIntensities(const std::vector<DoubleReal>& db_vec)))
{
    TEST_EXCEPTION(Exception::InvalidValue, test_mt.setSmoothedIntensities(smoothed_ints));
    smoothed_ints.pop_back();

    test_mt.setSmoothedIntensities(smoothed_ints);

    TEST_EQUAL(test_mt.getSmoothedIntensities().size(), smoothed_ints.size())
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

test_mt.setSmoothedIntensities(smoothed_ints);

START_SECTION((DoubleReal getIntensity(bool)))
{
    TEST_EXCEPTION(Exception::InvalidValue, test_mt.getIntensity(true));

    test_mt.estimateFWHM(true);

    DoubleReal smoothed_area = test_mt.getIntensity(true);
    TEST_REAL_SIMILAR(smoothed_area, 69460700);

    DoubleReal raw_area = test_mt.getIntensity(false);
    TEST_REAL_SIMILAR(raw_area, 69922872.7);
}
END_SECTION

/////

START_SECTION((DoubleReal getMaxIntensity(bool)))
{
    DoubleReal smoothed_maxint = test_mt.getMaxIntensity(true);
    TEST_REAL_SIMILAR(smoothed_maxint, 33000000.0);

    DoubleReal raw_maxint= test_mt.getMaxIntensity(false);
    TEST_REAL_SIMILAR(raw_maxint, 33329536.0);
}
END_SECTION

/////

START_SECTION((DoubleReal getMaxIntensity(bool) const ))
{
    const MassTrace test_mt_const(test_mt);
    DoubleReal smoothed_maxint = test_mt_const.getMaxIntensity(true);
    TEST_REAL_SIMILAR(smoothed_maxint, 33000000.0);

    DoubleReal raw_maxint= test_mt_const.getMaxIntensity(false);
    TEST_REAL_SIMILAR(raw_maxint, 33329536.0);
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


START_SECTION((DoubleReal getFWHM()))
{
    DoubleReal test_mt_fwhm = test_mt.getFWHM();

    TEST_REAL_SIMILAR(test_mt_fwhm, 4.01);
}
END_SECTION

/////

START_SECTION((DoubleReal getFWHM() const ))
{
    MassTrace test_mt_const(test_mt);

    DoubleReal test_mt_fwhm = test_mt_const.getFWHM();

    TEST_REAL_SIMILAR(test_mt_fwhm, 4.01);
}
END_SECTION

/////

START_SECTION((DoubleReal computeSmoothedPeakArea()))
{
    DoubleReal peak_area = test_mt.computeSmoothedPeakArea();

    TEST_REAL_SIMILAR(peak_area, 70129723.0)
}
END_SECTION

/////


START_SECTION((DoubleReal computePeakArea()))
{
    DoubleReal peak_area = test_mt.computePeakArea();

    TEST_REAL_SIMILAR(peak_area, 70164277.0)
}
END_SECTION

/////

START_SECTION((DoubleReal computePeakArea() const ))
{
    MassTrace test_mt_const(test_mt);
    DoubleReal peak_area = test_mt_const.computePeakArea();

    TEST_REAL_SIMILAR(peak_area, 70164277.0)
}
END_SECTION

/////

START_SECTION((DoubleReal computeFwhmAreaSmooth()))
{
    DoubleReal peak_area = test_mt.computeFwhmAreaSmooth();

    TEST_REAL_SIMILAR(peak_area, 69040000.0)
}
END_SECTION

/////

START_SECTION((DoubleReal computeFwhmArea()))
{
    DoubleReal peak_area = test_mt.computeFwhmArea();

    TEST_REAL_SIMILAR(peak_area, 69497153.0)
}
END_SECTION

/////

START_SECTION((DoubleReal computeFwhmAreaSmoothRobust()))
{
    DoubleReal peak_area = test_mt.computeFwhmAreaSmoothRobust();

    TEST_REAL_SIMILAR(peak_area, 69460700.0)
}
END_SECTION

/////

START_SECTION((DoubleReal computeFwhmAreaRobust()))
{
    DoubleReal peak_area = test_mt.computeFwhmAreaRobust();

    TEST_REAL_SIMILAR(peak_area, 69922872.67)
}
END_SECTION

/////

START_SECTION((Size findMaxByIntPeak(bool) const ))
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

START_SECTION((DoubleReal estimateFWHM(bool)))
{
    TEST_EXCEPTION(Exception::InvalidValue, test_mt2.estimateFWHM(true));
    TEST_EXCEPTION(Exception::InvalidValue, test_mt3.estimateFWHM(false));

    DoubleReal test_fwhm1 = test_mt.estimateFWHM(false);
    DoubleReal test_fwhm2 = test_mt.estimateFWHM(true);

    TEST_REAL_SIMILAR(test_fwhm1, 4.01);
    TEST_REAL_SIMILAR(test_fwhm2, 4.01);

    // TEST_EQUAL(test_mt.getFWHMScansNum(), 5);
}
END_SECTION

/////

START_SECTION((std::pair<Size, Size> getFWHMborders()))
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

std::vector<DoubleReal> double_smooth_ints(smoothed_ints);
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

    std::vector<DoubleReal> sm1(copy_mt.getSmoothedIntensities()), sm2(test_mt.getSmoothedIntensities());

    std::vector<DoubleReal>::const_iterator sm1_it = sm1.begin();

    for (std::vector<DoubleReal>::const_iterator sm2_it = sm2.begin(); sm2_it != sm2.end(); ++sm2_it)
    {
        TEST_EQUAL(*sm1_it, *sm2_it);
        ++sm1_it;
    }


    // TEST_EQUAL(copy_mt.getFWHMScansNum(), test_mt.getFWHMScansNum());
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


    // TEST_EQUAL(copy_mt.getFWHMScansNum(), test_mt.getFWHMScansNum());
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

    TEST_REAL_SIMILAR(test_mt.getCentroidRT(),155.22051);
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

    DoubleReal smooth_max_rt = test_mt.getCentroidRT();

    TEST_REAL_SIMILAR(smooth_max_rt, 155.24);
}
END_SECTION

/////

START_SECTION((void updateSmoothedWeightedMeanRT()))
{
    MassTrace raw_mt(peak_vec);

    TEST_EXCEPTION(Exception::InvalidValue, raw_mt.updateSmoothedWeightedMeanRT());

    test_mt.updateSmoothedWeightedMeanRT();

    DoubleReal smooth_max_rt = test_mt.getCentroidRT();

    TEST_REAL_SIMILAR(smooth_max_rt, 155.2389);
}
END_SECTION

/////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
