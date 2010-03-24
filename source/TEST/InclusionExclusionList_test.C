// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/InclusionExclusionList.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(InclusionExclusionList, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InclusionExclusionList* ptr = 0;
START_SECTION(InclusionExclusionList())
{
	ptr = new InclusionExclusionList();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~InclusionExclusionList())
{
	delete ptr;
}
END_SECTION

START_SECTION((void writeTargets(std::vector< FASTAFile::FASTAEntry > &fasta_entries, String &out_path, IntList &charges, String rt_model_path, DoubleReal rel_rt_window_size, bool rt_in_seconds, Size missed_cleavages=0)))
{
	// load data and write out file
	InclusionExclusionList list;
	vector<FASTAFile::FASTAEntry> entries;
	FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_1.fasta"),entries);
	IntList charges;
	charges<<2;
	String rt_model_path(OPENMS_GET_TEST_DATA_PATH("RTSimulation_absolut_rt.model"));
	DoubleReal rel_rt_window_size = 0.05;
	bool rt_in_seconds = true;
	Size missed_cleavages = 0;
	String out;
	NEW_TMP_FILE(out);
	// rt in seconds
	list.writeTargets(entries,out,charges,rt_model_path,rel_rt_window_size,rt_in_seconds,missed_cleavages);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_1_out.txt"),out)

	// rt in minutes
	rt_in_seconds = false;
	String out2;
	NEW_TMP_FILE(out2);
	list.writeTargets(entries,out2,charges,rt_model_path,rel_rt_window_size,rt_in_seconds,missed_cleavages);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_1_minutes_out.txt"),out2)
}
END_SECTION

START_SECTION((void writeTargets(FeatureMap<> &map, String &out_path, DoubleReal rel_rt_window_size, bool rt_in_seconds)))
{
  InclusionExclusionList list;
	FeatureMap<> map;
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_2.featureXML"),map);
	DoubleReal rel_rt_window_size = 0.05;
	bool rt_in_seconds = true;
	String out;
	NEW_TMP_FILE(out);
	list.writeTargets(map,out,rel_rt_window_size,rt_in_seconds);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_2_out.txt"),out)
	String out2;
	NEW_TMP_FILE(out2);
	rt_in_seconds = false;
	list.writeTargets(map,out2,rel_rt_window_size,rt_in_seconds);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_2_minutes_out.txt"),out2)
	
}
END_SECTION

START_SECTION((void writeTargets(std::vector< PeptideIdentification > &pep_ids, String &out_path, DoubleReal rel_rt_window_size, IntList &charges, bool rt_in_seconds)))
{
  InclusionExclusionList list;
	FeatureMap<> map;
	vector<PeptideIdentification> pep_ids;
	vector<ProteinIdentification> prot_ids;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3.IdXML"),prot_ids,pep_ids);
	DoubleReal rel_rt_window_size = 0.05;
	bool rt_in_seconds = true;
	IntList charges;
	charges<<2;
	String out;
	NEW_TMP_FILE(out);
	list.writeTargets(pep_ids,out,rel_rt_window_size,charges,rt_in_seconds);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3_out.txt"),out)
	String out2;
	NEW_TMP_FILE(out2);
	rt_in_seconds = false;
	list.writeTargets(pep_ids,out2,rel_rt_window_size,charges,rt_in_seconds);
	TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("InclusionExclusionList_3_minutes_out.txt"),out2)
  
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



