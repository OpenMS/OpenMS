// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricQuantifierStatistics, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricQuantifierStatistics* ptr = nullptr;
IsobaricQuantifierStatistics* null_ptr = nullptr;
START_SECTION(IsobaricQuantifierStatistics())
{
	ptr = new IsobaricQuantifierStatistics();
	TEST_NOT_EQUAL(ptr, null_ptr)
    
  TEST_EQUAL(ptr->channel_count, 0)
  TEST_EQUAL(ptr->iso_number_ms2_negative, 0)
  TEST_EQUAL(ptr->iso_number_reporter_negative, 0)
  TEST_EQUAL(ptr->iso_number_reporter_different, 0)
  TEST_EQUAL(ptr->iso_solution_different_intensity, 0.0)
  TEST_EQUAL(ptr->iso_total_intensity_negative, 0.0)
  TEST_EQUAL(ptr->number_ms2_total, 0)
  TEST_EQUAL(ptr->number_ms2_empty, 0)
  TEST_EQUAL(ptr->empty_channels.empty(), true)
}
END_SECTION

START_SECTION(~IsobaricQuantifierStatistics())
{
	delete ptr;
}
END_SECTION

START_SECTION((void reset()))
{
  IsobaricQuantifierStatistics stats;

  stats.channel_count = 4;
  stats.iso_number_ms2_negative = 10;
  stats.iso_number_reporter_negative = 20;
  stats.iso_number_reporter_different = 10;
  stats.iso_solution_different_intensity = 131.3;
  stats.iso_total_intensity_negative = 134.3;
  stats.number_ms2_total = 200;
  stats.number_ms2_empty = 3;
  stats.empty_channels[114] = 4;
  
  stats.reset();
  
  // check if reset worked properly
  TEST_EQUAL(stats.channel_count, 0)
  TEST_EQUAL(stats.iso_number_ms2_negative, 0)
  TEST_EQUAL(stats.iso_number_reporter_negative, 0)
  TEST_EQUAL(stats.iso_number_reporter_different, 0)
  TEST_EQUAL(stats.iso_solution_different_intensity, 0.0)
  TEST_EQUAL(stats.iso_total_intensity_negative, 0.0)
  TEST_EQUAL(stats.number_ms2_total, 0)
  TEST_EQUAL(stats.number_ms2_empty, 0)
  TEST_EQUAL(stats.empty_channels.empty(), true)
}
END_SECTION

START_SECTION((IsobaricQuantifierStatistics(const IsobaricQuantifierStatistics &other)))
{
  IsobaricQuantifierStatistics stats;

  stats.channel_count = 4;
  stats.iso_number_ms2_negative = 10;
  stats.iso_number_reporter_negative = 20;
  stats.iso_number_reporter_different = 10;
  stats.iso_solution_different_intensity = 131.3;
  stats.iso_total_intensity_negative = 134.3;
  stats.number_ms2_total = 200;
  stats.number_ms2_empty = 3;
  stats.empty_channels[114] = 4;
  
  IsobaricQuantifierStatistics stats2(stats);
  TEST_EQUAL(stats2.channel_count, 4)
  TEST_EQUAL(stats2.iso_number_ms2_negative, 10)
  TEST_EQUAL(stats2.iso_number_reporter_negative, 20)
  TEST_EQUAL(stats2.iso_number_reporter_different, 10)
  TEST_EQUAL(stats2.iso_solution_different_intensity, 131.3)
  TEST_EQUAL(stats2.iso_total_intensity_negative, 134.3)
  TEST_EQUAL(stats2.number_ms2_total, 200)
  TEST_EQUAL(stats2.number_ms2_empty, 3)
  TEST_EQUAL(stats2.empty_channels.find(114) != stats2.empty_channels.end(), true)
  TEST_EQUAL(stats2.empty_channels[114], 4)
}
END_SECTION

START_SECTION((IsobaricQuantifierStatistics& operator=(const IsobaricQuantifierStatistics &rhs)))
{
  IsobaricQuantifierStatistics stats;

  stats.channel_count = 4;
  stats.iso_number_ms2_negative = 10;
  stats.iso_number_reporter_negative = 20;
  stats.iso_number_reporter_different = 10;
  stats.iso_solution_different_intensity = 131.3;
  stats.iso_total_intensity_negative = 134.3;
  stats.number_ms2_total = 200;
  stats.number_ms2_empty = 3;
  stats.empty_channels[114] = 4;
  
  IsobaricQuantifierStatistics stats2;
  stats2 = stats;
  
  TEST_EQUAL(stats2.channel_count, 4)
  TEST_EQUAL(stats2.iso_number_ms2_negative, 10)
  TEST_EQUAL(stats2.iso_number_reporter_negative, 20)
  TEST_EQUAL(stats2.iso_number_reporter_different, 10)
  TEST_EQUAL(stats2.iso_solution_different_intensity, 131.3)
  TEST_EQUAL(stats2.iso_total_intensity_negative, 134.3)
  TEST_EQUAL(stats2.number_ms2_total, 200)
  TEST_EQUAL(stats2.number_ms2_empty, 3)
  TEST_EQUAL(stats2.empty_channels.find(114) != stats2.empty_channels.end(), true)
  TEST_EQUAL(stats2.empty_channels[114], 4)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
