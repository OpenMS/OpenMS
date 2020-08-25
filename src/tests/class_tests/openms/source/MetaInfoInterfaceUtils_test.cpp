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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/MetaInfoInterfaceUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MetaInfoInterfaceUtils, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((template<typename T_In, T_Out> static T_Out findCommonMetaKeys(const T::const_iterator& start, const T::const_iterator& end, const float min_frequency = 100.0)))
{
  vector<PeptideHit> hits; // some class derived from MetaInfoInterface
  for (Size i = 0; i < 10; ++i)
  {
    PeptideHit h;
    h.setMetaValue("commonMeta1", i);
    h.setMetaValue("commonMeta2", i);
    if (i % 2 == 0) 
    {
      h.setMetaValue("meta50pc", i);
    }
    hits.push_back(h);
  }
  hits.back().setMetaValue("metaSingle", "single");
  
  // common keys for ALL entries (i.e. 100% min_frequency)
  {
    std::vector<String> common = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::vector<String> >(hits.begin(), hits.end(), 100.0);
    TEST_EQUAL(common.size(), 2);
    ABORT_IF(common.size() != 2);
    TEST_EQUAL(common[0], "commonMeta1");
    TEST_EQUAL(common[1], "commonMeta2");
    
    // exceeds 100% --> should be corrected to 100% internally
    std::vector<String> common2 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::vector<String> >(hits.begin(), hits.end(), 1110.0);
    TEST_EQUAL(common==common2, true);
  }
  
  // occurrence of at least 50 (i.e. 50% min_frequency)
  {
    std::set<String> set50 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::set<String> >(hits.begin(), hits.end(), 50.0);
    TEST_EQUAL(set50.size(), 3);
    ABORT_IF(set50.size() != 3);
    std::set<String> set50_expected;
    set50_expected.insert("commonMeta1");
    set50_expected.insert("commonMeta2");
    set50_expected.insert("meta50pc");
    TEST_EQUAL(set50==set50_expected, true);
  }

  // ALL keys (i.e. 0% min_frequency)
  {
    std::set<String> set0 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::set<String> >(hits.begin(), hits.end(), 0.0);
    TEST_EQUAL(set0.size(), 4);
    ABORT_IF(set0.size() != 4);
    std::set<String> set0_expected;
    set0_expected.insert("commonMeta1");
    set0_expected.insert("commonMeta2");
    set0_expected.insert("meta50pc");
    set0_expected.insert("metaSingle");
    TEST_EQUAL(set0==set0_expected, true);

    // exceeds 0% --> should be corrected to 0% internally
    std::set<String> set0_2 = MetaInfoInterfaceUtils::findCommonMetaKeys<std::vector<PeptideHit>, std::set<String> >(hits.begin(), hits.end(), -10.0);
    TEST_EQUAL(set0==set0_2, true);
  }
  

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



