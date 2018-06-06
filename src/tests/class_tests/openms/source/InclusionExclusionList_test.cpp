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
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"

START_TEST(InclusionExclusionList, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InclusionExclusionList* ptr = nullptr;
InclusionExclusionList* nullPointer = nullptr;
START_SECTION(InclusionExclusionList())
{
	ptr = new InclusionExclusionList();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~InclusionExclusionList())
{
	delete ptr;
}
END_SECTION

START_SECTION((void writeTargets(const std::vector<FASTAFile::FASTAEntry>& fasta_entries, const String& out_path, const IntList& charges, const String rt_model_path)))
{
	// load data and write out file
	InclusionExclusionList list;
	vector<FASTAFile::FASTAEntry> entries;
	FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_1.fasta"),entries);
	IntList charges;
	charges.push_back(2);
	String rt_model_path(OPENMS_GET_TEST_DATA_PATH("RTSimulation_absolut_rt.model"));
  Param p = list.getParameters();
  p.setValue("missed_cleavages", 0);
  p.setValue("RT:unit", "seconds");
  list.setParameters(p);
	String out;
	NEW_TMP_FILE(out);
	// rt in seconds
	list.writeTargets(entries,out,charges,rt_model_path);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_1_out.txt"),out)

	// rt in minutes
	String out2;
	NEW_TMP_FILE(out2);
  p.setValue("RT:unit", "minutes");
  list.setParameters(p);
	list.writeTargets(entries,out2,charges,rt_model_path);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_1_minutes_out.txt"),out2)
}
END_SECTION

START_SECTION((void writeTargets(const FeatureMap& map, const String& out_path)))
{
  InclusionExclusionList list;
	FeatureMap map;
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_2.featureXML"),map);
  Param p = list.getParameters();
  p.setValue("missed_cleavages", 0);
  p.setValue("RT:unit", "seconds");
  list.setParameters(p);
	String out;
	NEW_TMP_FILE(out);
	list.writeTargets(map, out);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_2_out.txt"),out)
	String out2;
	NEW_TMP_FILE(out2);
  p.setValue("RT:unit", "minutes");
  list.setParameters(p);
	list.writeTargets(map, out2);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_2_minutes_out.txt"),out2)

  /// test clustering
  map.clear();
  Feature f;
  f.setCharge(1);
  f.setRT(100);

  // start putting data in...
  // close in m/z case
  f.setMZ(1000);
  map.push_back(f);
  f.setMZ(1000.00001);
  map.push_back(f);

  // non-overlapping RT case (singleton expected)
  f.setRT(150);
  map.push_back(f);

  // overlapping RT case
  f.setRT(1500);
  map.push_back(f);
  f.setRT(1510);
  map.push_back(f);

  // overlapping RT, but too far in m/z
  f.setRT(1520);
  f.setMZ(1001);
  map.push_back(f);

  p.setValue("merge:rt_tol", 0.0);
  p.setValue("merge:mz_tol", 10.0);
  p.setValue("merge:mz_tol_unit", "ppm");
  list.setParameters(p);
  list.writeTargets(map,out);
  TextFile tf;
  tf.load(out);

  TEST_EQUAL((tf.end() - tf.begin()), 4);

  // test exact m/z matching (no deviation allowed)
  {
  InclusionExclusionList list;
  Param p = list.getParameters();
  p.setValue("merge:rt_tol", 0.0);
  p.setValue("merge:mz_tol", 0.0);
  p.setValue("merge:mz_tol_unit", "ppm");
  list.setParameters(p);

  list.writeTargets(map,out);
  TextFile tf;
  tf.load(out);

  TEST_EQUAL((tf.end() - tf.begin()), 5);
  }

  // now test window overlap
  {
  InclusionExclusionList list;
  Param p = list.getParameters();
  p.setValue("merge:rt_tol", 11.0);
  p.setValue("merge:mz_tol", 0.0);
  p.setValue("merge:mz_tol_unit", "ppm");
  list.setParameters(p);
  list.writeTargets(map, out);
  TextFile tf;
  tf.load(out);

  TEST_EQUAL((tf.end() - tf.begin()), 5);
  }


}
END_SECTION

START_SECTION((void writeTargets(const std::vector<PeptideIdentification>& pep_ids, const String& out_path, const IntList& charges)))
{
  InclusionExclusionList list;
	FeatureMap map;
	vector<PeptideIdentification> pep_ids;
	vector<ProteinIdentification> prot_ids;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3.idXML"),prot_ids,pep_ids);
  Param p = list.getParameters();
  p.setValue("RT:unit", "seconds");
  list.setParameters(p);
	IntList charges;
	charges.push_back(2);
	String out;
	NEW_TMP_FILE(out);
	list.writeTargets(pep_ids,out,charges);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3_out.txt"),out)
	String out2;
	NEW_TMP_FILE(out2);
  p.setValue("RT:unit", "minutes");
  list.setParameters(p);
	list.writeTargets(pep_ids,out2,charges);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3_minutes_out.txt"),out2)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop

