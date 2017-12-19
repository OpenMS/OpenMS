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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
///////////////////////////

#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

using namespace OpenMS;

START_TEST(SwathWindowLoader, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SwathWindowLoader* ptr = nullptr;
SwathWindowLoader* nullPointer = nullptr;

START_SECTION(SwathWindowLoader())
  ptr = new SwathWindowLoader;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(virtual ~SwathWindowLoader())
    delete ptr;
END_SECTION

START_SECTION( static void readSwathWindows(const String & filename, std::vector<double> & swath_prec_lower_, std::vector<double> & swath_prec_upper_ ) )
{
  std::vector<double> swath_prec_lower;
  std::vector<double> swath_prec_upper;
  SwathWindowLoader::readSwathWindows(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"),
      swath_prec_lower, swath_prec_upper);

  TEST_EQUAL(swath_prec_lower.size(), swath_prec_upper.size())
  TEST_REAL_SIMILAR(swath_prec_lower[0], 400)
  TEST_REAL_SIMILAR(swath_prec_lower[1], 425)
  TEST_REAL_SIMILAR(swath_prec_lower[2], 450)
  TEST_REAL_SIMILAR(swath_prec_lower[3], 475)

  TEST_REAL_SIMILAR(swath_prec_upper[0], 425)
  TEST_REAL_SIMILAR(swath_prec_upper[1], 450)
  TEST_REAL_SIMILAR(swath_prec_upper[2], 475)
  TEST_REAL_SIMILAR(swath_prec_upper[3], 500)
}
END_SECTION

START_SECTION ( static void annotateSwathMapsFromFile(const std::string & filename, std::vector< OpenSwath::SwathMap >& swath_maps, bool doSort))
{
  // Initialize swath maps
  std::vector< OpenSwath::SwathMap > swath_maps(4);
  for (int i = 0; i < (int)swath_maps.size(); i++)
  {
    swath_maps[i].lower = 0;
    swath_maps[i].upper = 0;
    swath_maps[i].center = 0;
    swath_maps[i].ms1 = false;
  }

  SwathWindowLoader::annotateSwathMapsFromFile(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_maps, false);

  TEST_REAL_SIMILAR(swath_maps[0].lower, 400)
  TEST_REAL_SIMILAR(swath_maps[1].lower, 425)
  TEST_REAL_SIMILAR(swath_maps[2].lower, 450)
  TEST_REAL_SIMILAR(swath_maps[3].lower, 475)
                                 
  TEST_REAL_SIMILAR(swath_maps[0].upper, 425)
  TEST_REAL_SIMILAR(swath_maps[1].upper, 450)
  TEST_REAL_SIMILAR(swath_maps[2].upper, 475)
  TEST_REAL_SIMILAR(swath_maps[3].upper, 500)

  // Initialize swath maps
  std::vector< OpenSwath::SwathMap > swath_maps_sorted(4);
  for (int i = 0; i < (int)swath_maps_sorted.size(); i++)
  {
    swath_maps_sorted[i].lower = -i;
    swath_maps_sorted[i].upper = -i;
    swath_maps_sorted[i].center = -i;
    swath_maps_sorted[i].ms1 = false;
  }

  // test before
  TEST_REAL_SIMILAR(swath_maps_sorted[0].center, 0)
  TEST_REAL_SIMILAR(swath_maps_sorted[1].center, -1)
  TEST_REAL_SIMILAR(swath_maps_sorted[2].center, -2)
  TEST_REAL_SIMILAR(swath_maps_sorted[3].center, -3)

  SwathWindowLoader::annotateSwathMapsFromFile(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), 
      swath_maps_sorted, true);

  TEST_REAL_SIMILAR(swath_maps_sorted[0].lower, 400)
  TEST_REAL_SIMILAR(swath_maps_sorted[1].lower, 425)
  TEST_REAL_SIMILAR(swath_maps_sorted[2].lower, 450)
  TEST_REAL_SIMILAR(swath_maps_sorted[3].lower, 475)
                                        
  TEST_REAL_SIMILAR(swath_maps_sorted[0].upper, 425)
  TEST_REAL_SIMILAR(swath_maps_sorted[1].upper, 450)
  TEST_REAL_SIMILAR(swath_maps_sorted[2].upper, 475)
  TEST_REAL_SIMILAR(swath_maps_sorted[3].upper, 500)

  // should now be in reverse order
  TEST_REAL_SIMILAR(swath_maps_sorted[0].center, -3)
  TEST_REAL_SIMILAR(swath_maps_sorted[1].center, -2)
  TEST_REAL_SIMILAR(swath_maps_sorted[2].center, -1)
  TEST_REAL_SIMILAR(swath_maps_sorted[3].center, 0)

  // Test exceptions
  std::vector< OpenSwath::SwathMap > swath_maps_too_large(5);
  TEST_EXCEPTION(OpenMS::Exception::IllegalArgument, 
      SwathWindowLoader::annotateSwathMapsFromFile(
        OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_maps_too_large, true));

  std::vector< OpenSwath::SwathMap > swath_maps_too_small(3);
  TEST_EXCEPTION(OpenMS::Exception::IllegalArgument, 
      SwathWindowLoader::annotateSwathMapsFromFile(
        OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_maps_too_small, true));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
