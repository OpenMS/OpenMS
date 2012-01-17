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
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ElutionPeakDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ElutionPeakDetection* ptr = 0;
ElutionPeakDetection* null_ptr = 0;
START_SECTION(ElutionPeakDetection())
{
	ptr = new ElutionPeakDetection();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ElutionPeakDetection())
{
	delete ptr;
}
END_SECTION

MSExperiment<Peak1D> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ElutionPeakDetection_input1.mzML"), input);

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

MassTraceDetection test_mtd;
test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;

String mt_labels[] = {"T6", "T7", "T9", "T3", "T4", "T8", "T5", "T2", "T1"};
String split_labels[] = {"T6_0", "T6_1", "T6_2", "T6_3", "T7_0", "T7_1", "T7_2", "T7_3", "T7_4", "T9_0", "T9_1", "T9_2", "T3_0", "T3_1", "T3_2", "T3_3", "T3_4", "T3_5", "T3_6", "T3_7", "T3_8", "T3_9", "T3_10", "T3_11", "T3_12", "T4_0", "T4_1", "T4_2", "T8_0", "T8_1", "T5", "T2", "T1_0", "T1_1", "T1_2"};
String filt_labels[] = {"T6_0", "T6_1", "T7_3", "T3_1", "T3_2", "T3_3", "T3_5", "T3_7", "T3_8", "T3_9", "T4_1", "T1_0", "T6_3", "T3_0", "T3_4", "T3_11", "T4_0", "T4_2", "T8_1", "T6_2", "T7_1", "T7_2", "T9_0", "T3_10", "T8_0", "T7_0", "T5", "T9_1", "T2", "T1_1", "T1_2", "T3_6"};



START_SECTION((void detectPeaks(std::vector< MassTrace > &, std::vector< MassTrace > &)))
{
    TEST_EQUAL(output_mt.size(), 9);

    for (Size i = 0; i < output_mt.size(); ++i)
    {
        TEST_EQUAL(output_mt[i].getLabel(), mt_labels[i]);
    }


    test_epd.detectPeaks(output_mt, splitted_mt);

    // mass traces splitted to local peaks
    TEST_EQUAL(splitted_mt.size(), 35);

    for (Size i = 0; i < splitted_mt.size(); ++i)
    {
        TEST_EQUAL(splitted_mt[i].getLabel(), split_labels[i]);
    }
}
END_SECTION

START_SECTION((void detectPeaks(MassTrace &, std::vector< MassTrace > &)))
{
  NOT_TESTABLE; // see above
}
END_SECTION

START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
{
    test_epd.filterByPeakWidth(splitted_mt, filtered_mt);

    TEST_EQUAL(filtered_mt.size(), 32);

    for (Size i = 0; i < filtered_mt.size(); ++i)
    {
        TEST_EQUAL(filtered_mt[i].getLabel(), filt_labels[i]);
    }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



