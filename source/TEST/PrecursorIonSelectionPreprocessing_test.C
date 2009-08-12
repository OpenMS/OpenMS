// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PrecursorIonSelectionPreprocessing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PrecursorIonSelectionPreprocessing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PrecursorIonSelectionPreprocessing* ptr = 0;
START_SECTION(PrecursorIonSelectionPreprocessing())
	ptr = new PrecursorIonSelectionPreprocessing();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PrecursorIonSelectionPreprocessing())
	delete ptr;
END_SECTION

ptr = new PrecursorIonSelectionPreprocessing();

START_SECTION(PrecursorIonSelectionPreprocessing(const PrecursorIonSelectionPreprocessing &source))
	PrecursorIonSelectionPreprocessing copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PrecursorIonSelectionPreprocessing& operator=(const PrecursorIonSelectionPreprocessing &source))
  PrecursorIonSelectionPreprocessing copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

Param param;
param.setValue("precursor_mass_tolerance",0.9);
param.setValue("precursor_mass_tolerance_unit","Da");
param.setValue("missed_cleavages",0);
std::string tmp_filename;
NEW_TMP_FILE(tmp_filename);
param.setValue("preprocessing:preprocessed_db_path",tmp_filename);
ptr->setParameters(param);
ptr->dbPreprocessing(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_db.fasta"),true);
	
START_SECTION((const std::map<String,std::vector<DoubleReal> >& getProtMasses() const))
	std::map<String,std::vector<DoubleReal> > prot_map = ptr->getProtMasses();
	TEST_EQUAL(prot_map.size(), 3)
END_SECTION

START_SECTION((const std::vector<DoubleReal> & getMasses(String acc)))
	const std::vector<DoubleReal>& pep_masses= ptr->getMasses("P01008");
	TEST_EQUAL(pep_masses.size(), 14)
	TEST_REAL_SIMILAR(pep_masses[0],1356.68332791328)
	const std::vector<DoubleReal>& pep_masses2= ptr->getMasses("P02787");
  TEST_EQUAL(pep_masses2.size(), 19)
	TEST_REAL_SIMILAR(pep_masses2[0],306.159984588623)
END_SECTION

	
	START_SECTION(void dbPreprocessing(String db_path, bool save=true))
	std::map<String,std::vector<DoubleReal> > prot_map = ptr->getProtMasses();
	TEST_EQUAL(prot_map.size(), 3)
END_SECTION

START_SECTION(DoubleReal getWeight(DoubleReal mass))
  DoubleReal w = ptr->getWeight(147.113);
  TEST_REAL_SIMILAR(w,1)
END_SECTION

START_SECTION(void loadPreprocessing())
	PrecursorIonSelectionPreprocessing ldb;
  ldb.setParameters(param);
  ldb.loadPreprocessing();
  TEST_EQUAL(ldb.getProtMasses().size(),3)
  DoubleReal w = ldb.getWeight(147.113);
  TEST_REAL_SIMILAR(w,1)

	std::vector<DoubleReal> pep_masses_l = ldb.getMasses("P01008");
  std::vector<DoubleReal> pep_masses = ptr->getMasses("P01008");
  TEST_EQUAL(pep_masses_l.size(),pep_masses.size())
	TEST_REAL_SIMILAR(pep_masses_l[0],pep_masses[0])
END_SECTION

PrecursorIonSelectionPreprocessing rt_pt_pp;
rt_pt_pp.setParameters(param);
rt_pt_pp.dbPreprocessing(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_db.fasta"),
												 OPENMS_GET_TEST_DATA_PATH("RTSimulation.svm"),
												 OPENMS_GET_TEST_DATA_PATH("DetectabilitySimulation.svm"),false);
	
START_SECTION(void dbPreprocessing(String db_path,String rt_model_path,String dt_model_path,bool save=true))
  TEST_EQUAL(rt_pt_pp.getProtMasses().size(),3);
  DoubleReal w = rt_pt_pp.getWeight(147.113);
  TEST_REAL_SIMILAR(w,1)
	TEST_REAL_SIMILAR(rt_pt_pp.getRT("P01008",1),0.09798)
 	TEST_REAL_SIMILAR(rt_pt_pp.getPT("P01008",1),0.0402)	
END_SECTION	

START_SECTION(void DoubleReal getRT(String prot_id,Size peptide_index))
  TEST_REAL_SIMILAR(rt_pt_pp.getRT("P01008",1),0.09798)
END_SECTION
	
START_SECTION(void DoubleReal getPT(String prot_id,Size peptide_index))
  TEST_REAL_SIMILAR(rt_pt_pp.getPT("P01008",1),0.0402)
END_SECTION

START_SECTION(void DoubleReal getRTWeight(String prot_id, Size peptide_index,DoubleReal meas_rt))
	TEST_REAL_SIMILAR(rt_pt_pp.getRTWeight("P01008",1,0.2),99.999)
END_SECTION

	
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



