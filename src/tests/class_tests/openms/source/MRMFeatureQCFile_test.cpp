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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

#include <OpenMS/FORMAT/MRMFeatureQCFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeatureQCFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureQCFile* ptr = nullptr;
MRMFeatureQCFile* nullPointer = nullptr;

START_SECTION(MRMFeatureQCFile())
{
  ptr = new MRMFeatureQCFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureQCFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, MRMFeatureQC& mrmfqc, const bool is_component_group) const)
{
  MRMFeatureQCFile mrmfqcfile;
  MRMFeatureQC mrmfqc;
  mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_1.csv"), mrmfqc, false); // components file
  mrmfqcfile.load(OPENMS_GET_TEST_DATA_PATH("MRMFeatureQCFile_2.csv"), mrmfqc, true); // component groups file
  TEST_EQUAL(mrmfqc.component_qcs[0].component_name, "component1");
  TEST_EQUAL(mrmfqc.component_qcs[1].component_name, "component2");
  TEST_EQUAL(mrmfqc.component_qcs[2].component_name, "component3");
  TEST_REAL_SIMILAR(mrmfqc.component_qcs[0].meta_value_qc["sn_score"].second, 10.0);
  TEST_REAL_SIMILAR(mrmfqc.component_qcs[1].meta_value_qc["sn_score"].second, 20.0);
  TEST_REAL_SIMILAR(mrmfqc.component_qcs[2].meta_value_qc["sn_score"].second, 50.0);
  TEST_EQUAL(mrmfqc.component_group_qcs[0].component_group_name, "componentGroup1");
  TEST_EQUAL(mrmfqc.component_group_qcs[1].component_group_name, "componentGroup2");
  TEST_EQUAL(mrmfqc.component_group_qcs[2].component_group_name, "componentGroup3");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

