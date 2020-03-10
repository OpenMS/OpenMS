// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConversionHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((template < typename FeatureT > static void convert(UInt64 const input_map_index, FeatureMap< FeatureT > const &input_map, ConsensusMap &output_map, Size n=-1)))
{

  FeatureMap fm;
  Feature f;
  for ( UInt i = 0; i < 3; ++i )
  {
    f.setRT(i*77.7);
    f.setMZ(i+100.35);
    f.setUniqueId(i*33+17);
    fm.push_back(f);
  }
  ConsensusMap cm;
  MapConversion::convert(33,fm,cm);

  TEST_EQUAL(cm.size(),3);
  TEST_EQUAL(cm.getColumnHeaders()[33].size,3);
  for ( UInt i = 0; i < 3; ++i )
  {
    TEST_EQUAL(cm[i].size(),1);
    TEST_EQUAL(cm[i].begin()->getMapIndex(),33);
    TEST_EQUAL(cm[i].begin()->getUniqueId(),i*33+17);
    TEST_REAL_SIMILAR(cm[i].begin()->getRT(),i*77.7);
    TEST_REAL_SIMILAR(cm[i].begin()->getMZ(),i+100.35);
  }

cm.clear();
MapConversion::convert(33,fm,cm,2);
TEST_EQUAL(cm.size(),2);
TEST_EQUAL(cm.getColumnHeaders()[33].size,3);

}
END_SECTION

/////

// Prepare data
PeakMap mse;
{
  MSSpectrum mss;
  Peak1D p;
  for ( UInt m = 0; m < 3; ++m )
  {
    mss.clear(true);
    for ( UInt i = 0; i < 4; ++i )
    {
      p.setMZ( 10* m + i + 100.35);
      p.setIntensity( 900 + 7*m + 5*i );
      mss.push_back(p);
    }
    mse.addSpectrum(mss);
    mse.getSpectra().back().setRT(m*5);
  }
}

START_SECTION((static void convert(UInt64 const input_map_index, PeakMap & input_map, ConsensusMap& output_map, Size n = -1)))
{

  ConsensusMap cm;

  MapConversion::convert(33,mse,cm,8);

  TEST_EQUAL(cm.size(),8);

  for ( UInt i = 0; i < cm.size(); ++i)
  {
    STATUS("\n" << i << ": " << cm[i] );
  }

  TEST_EQUAL(cm.back().getIntensity(),912);

}
END_SECTION

/////

ConsensusMap cm;
MapConversion::convert(33,mse,cm,8);

START_SECTION((template < typename FeatureT > static void convert(ConsensusMap const &input_map, const bool keep_uids, FeatureMap< FeatureT > &output_map)))
{
    FeatureMap out_fm;
    MapConversion::convert(cm, true, out_fm);

    TEST_EQUAL(cm.getUniqueId(), out_fm.getUniqueId());
    TEST_EQUAL(cm.getProteinIdentifications().size(), out_fm.getProteinIdentifications().size());
    TEST_EQUAL(cm.getUnassignedPeptideIdentifications().size(), out_fm.getUnassignedPeptideIdentifications().size());
    TEST_EQUAL(cm.size(), out_fm.size());

    for (Size i = 0; i < cm.size(); ++i)
    {
        TEST_EQUAL(cm[i], out_fm[i]);
    }

    out_fm.clear();
    MapConversion::convert(cm, false, out_fm);
    TEST_NOT_EQUAL(cm.getUniqueId(), out_fm.getUniqueId());

    for (Size i = 0; i < cm.size(); ++i)
    {
        TEST_REAL_SIMILAR(cm[i].getRT(), out_fm[i].getRT());
        TEST_REAL_SIMILAR(cm[i].getMZ(), out_fm[i].getMZ());
        TEST_REAL_SIMILAR(cm[i].getIntensity(), out_fm[i].getIntensity());

        TEST_NOT_EQUAL(cm[i].getUniqueId(), out_fm[i].getUniqueId());
    }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



