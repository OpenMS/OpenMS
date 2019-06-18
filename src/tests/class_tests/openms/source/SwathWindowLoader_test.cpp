// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
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

START_SECTION( static void readSwathWindows(const std::string& filename, std::vector<double>& swath_prec_lower, std::vector<double>& swath_prec_upper) )
{
  std::vector<double> swath_prec_lower;
  std::vector<double> swath_prec_upper;
  SwathWindowLoader::readSwathWindows(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_prec_lower, swath_prec_upper);

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

START_SECTION ( static void annotateSwathMapsFromFile(const std::string& filename, std::vector<OpenSwath::SwathMap>& swath_maps, bool do_sort, bool force) )
{
  // pretend this is given in the raw data:
  const std::vector< OpenSwath::SwathMap > swath_maps = {
   {399, 426, 1, false},
   {424, 451, 2, false},
   {450, 475, 3, false}, // matches exacly (no overlap), but will be ok
   {474, 501, 4, false}
  };

  // copy to feed into function
  std::vector< OpenSwath::SwathMap > swath_maps_test = swath_maps;


  SwathWindowLoader::annotateSwathMapsFromFile(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_maps_test, false, false);

  TEST_REAL_SIMILAR(swath_maps_test[0].lower, 400)
  TEST_REAL_SIMILAR(swath_maps_test[1].lower, 425)
  TEST_REAL_SIMILAR(swath_maps_test[2].lower, 450)
  TEST_REAL_SIMILAR(swath_maps_test[3].lower, 475)
                                 
  TEST_REAL_SIMILAR(swath_maps_test[0].upper, 425)
  TEST_REAL_SIMILAR(swath_maps_test[1].upper, 450)
  TEST_REAL_SIMILAR(swath_maps_test[2].upper, 475)
  TEST_REAL_SIMILAR(swath_maps_test[3].upper, 500)

  ///////////////
  // test sorting (start inverted)
  std::vector< OpenSwath::SwathMap > swath_maps_inv = swath_maps;
  // invert
  std::sort(swath_maps_inv.rbegin(), swath_maps_inv.rend(), [](const OpenSwath::SwathMap& a, const OpenSwath::SwathMap& b) { return a.lower < b.lower; });

  swath_maps_test = swath_maps_inv;
  // test before
  TEST_EQUAL(swath_maps_test[0].lower, 474)
  TEST_EQUAL(swath_maps_test[0].center, 4)
  TEST_EQUAL(swath_maps_test[3].lower, 399)
  TEST_EQUAL(swath_maps_test[3].center, 1)

  SwathWindowLoader::annotateSwathMapsFromFile(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), 
      swath_maps_test, true, false);

  TEST_REAL_SIMILAR(swath_maps_test[0].lower, 400)
  TEST_REAL_SIMILAR(swath_maps_test[1].lower, 425)
  TEST_REAL_SIMILAR(swath_maps_test[2].lower, 450)
  TEST_REAL_SIMILAR(swath_maps_test[3].lower, 475)
                                        
  TEST_REAL_SIMILAR(swath_maps_test[0].upper, 425)
  TEST_REAL_SIMILAR(swath_maps_test[1].upper, 450)
  TEST_REAL_SIMILAR(swath_maps_test[2].upper, 475)
  TEST_REAL_SIMILAR(swath_maps_test[3].upper, 500)

  // should now be in original order
  TEST_REAL_SIMILAR(swath_maps_test[0].center, 1)
  TEST_REAL_SIMILAR(swath_maps_test[1].center, 2)
  TEST_REAL_SIMILAR(swath_maps_test[2].center, 3)
  TEST_REAL_SIMILAR(swath_maps_test[3].center, 4)

  ///////////////////////////////////
  // Test exceptions
  std::vector< OpenSwath::SwathMap > swath_maps_too_large(5);
  TEST_EXCEPTION(OpenMS::Exception::IllegalArgument,
      SwathWindowLoader::annotateSwathMapsFromFile(
        OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_maps_too_large, true, false));

  std::vector< OpenSwath::SwathMap > swath_maps_too_small(3);
  TEST_EXCEPTION(OpenMS::Exception::IllegalArgument,
      SwathWindowLoader::annotateSwathMapsFromFile(
        OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), swath_maps_too_small, true, false));

  // wrong order && no sorting --> fail
  swath_maps_test = swath_maps_inv;
  // test before
  TEST_EQUAL(swath_maps_test[0].lower, 474)
  TEST_EQUAL(swath_maps_test[0].center, 4)
  TEST_EQUAL(swath_maps_test[3].lower, 399)
  TEST_EQUAL(swath_maps_test[3].center, 1)
  TEST_EXCEPTION(OpenMS::Exception::IllegalArgument,
      SwathWindowLoader::annotateSwathMapsFromFile(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), 
                                                   swath_maps_test, false, false));
  // wrong order && no sorting && force --> ok
  swath_maps_test = swath_maps_inv;
  SwathWindowLoader::annotateSwathMapsFromFile(OPENMS_GET_TEST_DATA_PATH("SwathWindowFile.txt"), 
                                               swath_maps_test, false, true);
  // overwritten windows (not narrowing the range since in wrong order )... but ok, due to -force
  TEST_REAL_SIMILAR(swath_maps_test[0].lower, 400)
  TEST_REAL_SIMILAR(swath_maps_test[1].lower, 425)
  TEST_REAL_SIMILAR(swath_maps_test[2].lower, 450)
  TEST_REAL_SIMILAR(swath_maps_test[3].lower, 475)
                                        
  TEST_REAL_SIMILAR(swath_maps_test[0].upper, 425)
  TEST_REAL_SIMILAR(swath_maps_test[1].upper, 450)
  TEST_REAL_SIMILAR(swath_maps_test[2].upper, 475)
  TEST_REAL_SIMILAR(swath_maps_test[3].upper, 500)

  // should still be in reverse order
  TEST_REAL_SIMILAR(swath_maps_test[0].center, 4)
  TEST_REAL_SIMILAR(swath_maps_test[1].center, 3)
  TEST_REAL_SIMILAR(swath_maps_test[2].center, 2)
  TEST_REAL_SIMILAR(swath_maps_test[3].center, 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
