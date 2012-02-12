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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassTraceDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassTraceDetection* ptr = 0;
MassTraceDetection* null_ptr = 0;
START_SECTION(MassTraceDetection())
{
    ptr = new MassTraceDetection();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MassTraceDetection())
{
    delete ptr;
}
END_SECTION

MassTraceDetection test_mtd;

START_SECTION((void updateIterativeWeightedMeanMZ(const DoubleReal &, const DoubleReal &, DoubleReal &, DoubleReal &, DoubleReal &)))
{
    DoubleReal centroid_mz(150.22), centroid_int(25000000);
    DoubleReal new_mz1(150.34), new_int1(23043030);
    DoubleReal new_mz2(150.11), new_int2(1932392);

    std::vector<DoubleReal> mzs, ints;
    mzs.push_back(centroid_mz);
    mzs.push_back(new_mz1);
    mzs.push_back(new_mz2);
    ints.push_back(centroid_int);
    ints.push_back(new_int1);
    ints.push_back(new_int2);

    DoubleReal total_weight1(centroid_int + new_int1);
    DoubleReal total_weight2(centroid_int + new_int1 + new_int2);

    DoubleReal wmean1((centroid_mz * centroid_int + new_mz1 * new_int1)/total_weight1);
    DoubleReal wmean2((centroid_mz * centroid_int + new_mz1 * new_int1 + new_mz2 * new_int2)/total_weight2);

    DoubleReal prev_count(centroid_mz * centroid_int);
    DoubleReal prev_denom(centroid_int);

    test_mtd.updateIterativeWeightedMeanMZ(new_mz1, new_int1, centroid_mz, prev_count, prev_denom);

    TEST_REAL_SIMILAR(centroid_mz, wmean1);

    test_mtd.updateIterativeWeightedMeanMZ(new_mz2, new_int2, centroid_mz, prev_count, prev_denom);

    TEST_REAL_SIMILAR(centroid_mz, wmean2);

}
END_SECTION

// load a mzML file for testing the algorithm
MSExperiment<Peak1D> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_input1.mzML"),input);

Size exp_mt_lengths[3] = {16, 31, 93};
DoubleReal exp_mt_rts[3] = {346.836, 346.881, 347.883};
DoubleReal exp_mt_mzs[3] = {439.27594, 438.27241, 437.26675};
DoubleReal exp_mt_ints[3] = {474.28, 2560.47, 12788.44};

std::vector<MassTrace> output_mt;

START_SECTION((void run(const MSExperiment< Peak1D > &, std::vector< MassTrace > &)))
{
    test_mtd.run(input, output_mt);

    TEST_EQUAL(output_mt.size(), 3);

    for (Size i = 0; i < output_mt.size(); ++i)
    {
        TEST_EQUAL(output_mt[i].getSize(), exp_mt_lengths[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidRT(), exp_mt_rts[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidMZ(), exp_mt_mzs[i]);
        TEST_REAL_SIMILAR(output_mt[i].computePeakArea(), exp_mt_ints[i]);
    }
}
END_SECTION

std::vector<MassTrace> filt;

START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
{
    test_mtd.filterByPeakWidth(output_mt, filt);

    TEST_EQUAL(output_mt.size(), filt.size());

    for (Size i = 0; i < output_mt.size(); ++i)
    {
        TEST_EQUAL(output_mt[i].getFWHMScansNum(), filt[i].getFWHMScansNum());
    }
}
END_SECTION

MSExperiment<Peak1D>::ConstAreaIterator mt_it1 = input.areaBeginConst(342.0, 375.0, 437.1, 437.4);
MSExperiment<Peak1D>::ConstAreaIterator mt_it2 = input.areaBeginConst(342.0, 375.0, 438.2, 438.4);
MSExperiment<Peak1D>::ConstAreaIterator mt_it3 = input.areaBeginConst(342.0, 375.0, 439.2, 439.4);

std::vector<MassTrace> found_mtraces;

MSExperiment<Peak1D>::ConstAreaIterator mt_end = input.areaEndConst();

START_SECTION((void run(MSExperiment< Peak1D >::ConstAreaIterator &begin, MSExperiment< Peak1D >::ConstAreaIterator &end, std::vector< MassTrace > &found_masstraces)))
{

    test_mtd.run(mt_it1, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[2]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[2]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[2]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[2]);

    found_mtraces.clear();


    test_mtd.run(mt_it2, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[1]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[1]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[1]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[1]);

    found_mtraces.clear();


    test_mtd.run(mt_it3, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[0]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[0]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[0]);

    found_mtraces.clear();
}
END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



