// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(InclusionExclusionList, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InclusionExclusionList* ptr = 0;
InclusionExclusionList* nullPointer = 0;
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
	charges<<2;
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

START_SECTION((void writeTargets(const FeatureMap<>& map, const String& out_path)))
{
  InclusionExclusionList list;
	FeatureMap<> map;
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

  TEST_EQUAL(tf.size(), 4);
  for (Size ii=0; ii<tf.size(); ++ii)
  {

    std::cout << tf[ii] << "\n";
  }

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

  TEST_EQUAL(tf.size(), 5);
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

  TEST_EQUAL(tf.size(), 5);
  }


}
END_SECTION

START_SECTION((void writeTargets(const std::vector<PeptideIdentification>& pep_ids, const String& out_path, const IntList& charges)))
{
  InclusionExclusionList list;
	FeatureMap<> map;
	vector<PeptideIdentification> pep_ids;
	vector<ProteinIdentification> prot_ids;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3.IdXML"),prot_ids,pep_ids);
	DoubleReal rel_rt_window_size = 0.05;
  Param p = list.getParameters();
  p.setValue("RT:unit", "seconds");
  list.setParameters(p);
	IntList charges;
	charges<<2;
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



