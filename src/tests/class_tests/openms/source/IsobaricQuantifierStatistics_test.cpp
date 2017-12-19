// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
