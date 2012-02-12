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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFindingMetabo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFindingMetabo* ptr = 0;
FeatureFindingMetabo* null_ptr = 0;
START_SECTION(FeatureFindingMetabo())
{
	ptr = new FeatureFindingMetabo();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFindingMetabo())
{
	delete ptr;
}
END_SECTION

// load a mzML file for testing the algorithm
MSExperiment<Peak1D> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_input1.mzML"), input);

FeatureMap<> exp_fm, test_fm;
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_output1.featureXML"), exp_fm);
exp_fm.sortByMZ();

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

MassTraceDetection test_mtd;
test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;
test_epd.detectPeaks(output_mt, splitted_mt);
test_epd.filterByPeakWidth(splitted_mt, filtered_mt);


START_SECTION((void run(std::vector< MassTrace > &, FeatureMap<> &)))
{
    FeatureFindingMetabo test_ffm;
    test_ffm.run(filtered_mt, test_fm);
    test_fm.sortByMZ();

    TEST_EQUAL(exp_fm.size(), test_fm.size());

    for (Size i = 0; i < exp_fm.size(); ++i)
    {
        TEST_EQUAL(exp_fm[i].getMetaValue(3), test_fm[i].getMetaValue(3));
        TEST_REAL_SIMILAR(exp_fm[i].getRT(), test_fm[i].getRT());
        TEST_REAL_SIMILAR(exp_fm[i].getMZ(), test_fm[i].getMZ());
        TEST_REAL_SIMILAR(exp_fm[i].getIntensity(), test_fm[i].getIntensity());
    }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



