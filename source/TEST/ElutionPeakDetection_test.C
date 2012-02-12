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
String split_labels[] = {"T6.1", "T6.2", "T6.3", "T6.4", "T7.1", "T7.2", "T7.3", "T7.4", "T7.5", "T9.1", "T9.2", "T9.3", "T3.1", "T3.2", "T3.3", "T3.4", "T3.5", "T3.6", "T3.7", "T3.8", "T3.9", "T3.10", "T3.11", "T3.12", "T3.13", "T4.1", "T4.2", "T4.3", "T8.1", "T8.2", "T5", "T2", "T1.1", "T1.2", "T1.3"};
String filt_labels[] = {"T6.1", "T6.2", "T7.4", "T3.2", "T3.3", "T3.4", "T3.6", "T3.8", "T3.9", "T3.10", "T4.2", "T1.1", "T6.4", "T3.1", "T3.5", "T3.12", "T4.1", "T4.3", "T8.2", "T6.3", "T7.2", "T7.3", "T9.1", "T3.11", "T8.1", "T7.1", "T5", "T9.2", "T2", "T1.2", "T1.3", "T3.7"};


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



