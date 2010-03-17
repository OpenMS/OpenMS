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
// $Maintainer: Alexandra Zerck$
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PrecursorIonSelection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PrecursorIonSelection* ptr = 0;
START_SECTION(PrecursorIonSelection())
  ptr = new PrecursorIonSelection();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PrecursorIonSelection())
  delete ptr;
END_SECTION

ptr = new PrecursorIonSelection();

START_SECTION(PrecursorIonSelection(const PrecursorIonSelection& source))
  ptr->setMaxScore(23.5);
  PrecursorIonSelection copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
  TEST_REAL_SIMILAR(copy.getMaxScore(), ptr->getMaxScore())
END_SECTION

START_SECTION(const DoubleReal& getMaxScore() const)
 TEST_REAL_SIMILAR(ptr->getMaxScore(),23.5)
END_SECTION  

START_SECTION(void setMaxScore(const DoubleReal& max_score))
  ptr->setMaxScore(24.5);
  TEST_REAL_SIMILAR(ptr->getMaxScore(),24.5)
END_SECTION

std::vector<ProteinIdentification> prot_ids;
std::vector<PeptideIdentification> pep_ids;
String document_id;
IdXMLFile file;
file.load(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelection_ids.IdXML"),prot_ids,pep_ids, document_id);

FeatureMap<> features,next_features;
FeatureXMLFile f_file;
f_file.load(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelection_features.featureXML"),features);
START_SECTION(void sortByTotalScore(FeatureMap<>& features))
  ptr->sortByTotalScore(features);
  TEST_REAL_SIMILAR((DoubleReal)features[0].getMetaValue("msms_score"),49485.75)
END_SECTION

START_SECTION(void getNextPrecursors(FeatureMap<>& features,FeatureMap<>& next_features,UInt number))
  ptr->getNextPrecursors(features,next_features,2);
  TEST_EQUAL(next_features.size(),2)
  TEST_REAL_SIMILAR((DoubleReal)next_features[0].getMetaValue("msms_score"),49485.75)
  TEST_REAL_SIMILAR((DoubleReal)next_features[1].getMetaValue("msms_score"),47365)
END_SECTION

PrecursorIonSelectionPreprocessing preprocessing;
Param param;
param.setValue("precursor_mass_tolerance",0.05);
param.setValue("precursor_mass_tolerance_unit","Da");
param.setValue("missed_cleavages",1);
param.setValue("preprocessing:preprocessed_db_path",OPENMS_GET_TEST_DATA_PATH(""));
preprocessing.setParameters(param);
preprocessing.dbPreprocessing(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelection_db.fasta"),false);

param.setValue("max_iteration",10);
param.setValue("type","IPS");
param.setValue("min_pep_ids",2);
param.remove("preprocessing:preprocessed_db_path");
ptr->setParameters(param);
next_features.clear(true);

START_SECTION(void rescore(FeatureMap<>& features,std::vector<PeptideIdentification>& new_pep_ids,std::vector<ProteinIdentification>& prot_ids,PrecursorIonSelectionPreprocessing& preprocessed_db, bool check_meta_values=true))
  ptr->rescore(features,pep_ids,prot_ids,preprocessing,false);
  ptr->getNextPrecursors(features,next_features,1);
  TEST_REAL_SIMILAR(next_features[0].getMetaValue("msms_score"),46365.5)
END_SECTION

START_SECTION( void simulateRun(FeatureMap<>& features,std::vector<PeptideIdentification>& pep_ids,std::vector<ProteinIdentification>& prot_ids,PrecursorIonSelectionPreprocessing& preprocessed_db,UInt step_size, String path))
  ptr->reset();
	features.clear(true);
  f_file.load(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelection_features.featureXML"),features);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  ptr->simulateRun(features,pep_ids,prot_ids,preprocessing,1,tmp_filename);
  ptr->sortByTotalScore(features);
  TEST_EQUAL(features[20].getMetaValue("shifted"),"both")
  TEST_REAL_SIMILAR(features[20].getMetaValue("msms_score"),27574.40625)
END_SECTION

START_SECTION((const std::map<String,std::set<String> >& getPeptideProteinCounter()))
	 TEST_EQUAL(ptr->getPeptideProteinCounter().size(),1)
END_SECTION
	
START_SECTION((void reset()))
	ptr->reset();
  TEST_EQUAL(ptr->getPeptideProteinCounter().size(),0)
END_SECTION

	
	
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



