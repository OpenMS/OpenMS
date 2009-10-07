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
// $Maintainer: Marc Sturm, Chris Bielow, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

DRange<1> makeRange(DoubleReal a, DoubleReal b)
{
  DPosition<1> pa(a), pb(b);
  return DRange<1>(pa, pb);
}

///////////////////////////

START_TEST(FeatureXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureXMLFile* ptr = 0;
START_SECTION((FeatureXMLFile()))
	ptr = new FeatureXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~FeatureXMLFile()))
	delete ptr;
END_SECTION

START_SECTION((void load(String filename, FeatureMap<>& feature_map)))
	TOLERANCE_ABSOLUTE(0.01)

	FeatureMap<> e;
	FeatureXMLFile dfmap_file;

	//test exception
	TEST_EXCEPTION( Exception::FileNotFound , dfmap_file.load("dummy/dummy.MzData",e) )

	// real test
	dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),e);
	TEST_EQUAL(e.getIdentifier(),"lsid");

	//test DocumentIdentifier addition
	TEST_STRING_EQUAL(e.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"));
	TEST_STRING_EQUAL(FileHandler::typeToName(e.getLoadedFileType()),"FeatureXML");

	TEST_EQUAL(e.size(),2)
	TEST_REAL_SIMILAR(e[0].getRT(), 25)
	TEST_REAL_SIMILAR(e[0].getMZ(), 0)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 300)
	TEST_EQUAL(e[0].getMetaValue("stringparametername"),"stringparametervalue")
	TEST_EQUAL((UInt)e[0].getMetaValue("intparametername"),4)
	TEST_REAL_SIMILAR((DoubleReal)e[0].getMetaValue("floatparametername"),4.551)
	TEST_REAL_SIMILAR(e[1].getRT(), 0)
	TEST_REAL_SIMILAR(e[1].getMZ(), 35)
	TEST_REAL_SIMILAR(e[1].getIntensity(), 500)
	//data processing
	TEST_EQUAL(e.getDataProcessing().size(),2)
	TEST_STRING_EQUAL(e.getDataProcessing()[0].getSoftware().getName(),"Software1")
	TEST_STRING_EQUAL(e.getDataProcessing()[0].getSoftware().getVersion(),"0.91a")
	TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().size(),1)
	TEST_EQUAL(e.getDataProcessing()[0].getProcessingActions().count(DataProcessing::DEISOTOPING),1)
	TEST_STRING_EQUAL(e.getDataProcessing()[0].getMetaValue("name"),"dataProcessing")
	TEST_STRING_EQUAL(e.getDataProcessing()[1].getSoftware().getName(),"Software2")
	TEST_STRING_EQUAL(e.getDataProcessing()[1].getSoftware().getVersion(),"0.92a")
	TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().size(),2)
	TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::SMOOTHING),1)
	TEST_EQUAL(e.getDataProcessing()[1].getProcessingActions().count(DataProcessing::BASELINE_REDUCTION),1)
	//protein identifications
	TEST_EQUAL(e.getProteinIdentifications().size(),2)
	TEST_EQUAL(e.getProteinIdentifications()[0].getHits().size(),2)
	TEST_EQUAL(e.getProteinIdentifications()[0].getHits()[0].getSequence(),"ABCDEFG")
	TEST_EQUAL(e.getProteinIdentifications()[0].getHits()[1].getSequence(),"HIJKLMN")
	TEST_EQUAL(e.getProteinIdentifications()[1].getHits().size(),1)
	TEST_EQUAL(e.getProteinIdentifications()[1].getHits()[0].getSequence(),"OPQREST")
	//peptide identifications
	TEST_EQUAL(e[0].getPeptideIdentifications().size(),2)
	TEST_EQUAL(e[0].getPeptideIdentifications()[0].getHits().size(),1)
	TEST_EQUAL(e[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),"A")
	TEST_EQUAL(e[0].getPeptideIdentifications()[1].getHits().size(),2)
	TEST_EQUAL(e[0].getPeptideIdentifications()[1].getHits()[0].getSequence(),"C")
	TEST_EQUAL(e[0].getPeptideIdentifications()[1].getHits()[1].getSequence(),"D")
	TEST_EQUAL(e[1].getPeptideIdentifications().size(),1)
	TEST_EQUAL(e[1].getPeptideIdentifications()[0].getHits().size(),1)
	TEST_EQUAL(e[1].getPeptideIdentifications()[0].getHits()[0].getSequence(),"E")
	//unassigned peptide identifications
	TEST_EQUAL(e.getUnassignedPeptideIdentifications().size(),2)
	TEST_EQUAL(e.getUnassignedPeptideIdentifications()[0].getHits().size(),1)
	TEST_EQUAL(e.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(),"F")
	TEST_EQUAL(e.getUnassignedPeptideIdentifications()[1].getHits().size(),2)
	TEST_EQUAL(e.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(),"G")
	TEST_EQUAL(e.getUnassignedPeptideIdentifications()[1].getHits()[1].getSequence(),"H")
	
	//test if loading a second file works (initialization)
	FeatureMap<> e2;
	dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),e2);
	TEST_EQUAL(e==e2,true)

	//test of old file with mzData description (version 1.2)
	//here only the downward-compatibility of the new parser is tested
	//no exception should be thrown
	dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_3_old.featureXML"),e);
	TEST_EQUAL(e.size(),1)

	//PeakFileOptions tests
	dfmap_file.getOptions().setRTRange(makeRange(0, 10));
	dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_SIMILAR(e[0].getRT(), 0)
	TEST_REAL_SIMILAR(e[0].getMZ(), 35)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 500)

	dfmap_file.getOptions() = PeakFileOptions();
	dfmap_file.getOptions().setMZRange(makeRange(10, 50));
	dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_SIMILAR(e[0].getRT(), 0)
	TEST_REAL_SIMILAR(e[0].getMZ(), 35)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 500)

	dfmap_file.getOptions() = PeakFileOptions();
	dfmap_file.getOptions().setIntensityRange(makeRange(400, 600));
	dfmap_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),e);
	TEST_EQUAL(e.size(),1)
	TEST_REAL_SIMILAR(e[0].getRT(), 0)
	TEST_REAL_SIMILAR(e[0].getMZ(), 35)
	TEST_REAL_SIMILAR(e[0].getIntensity(), 500)


END_SECTION

START_SECTION((void store(String filename, const FeatureMap<> &feature_map)))
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  FeatureMap<> map, map2;
  FeatureXMLFile f;

  f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),map);
  f.store(tmp_filename, map);
  f.load(tmp_filename, map2);
  TEST_EQUAL(map==map2, true)
END_SECTION

START_SECTION((PeakFileOptions& getOptions()))
	FeatureXMLFile f;
  FeatureMap<> e;
	f.getOptions().setRTRange(makeRange(1.5, 4.5));
	f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"),e);
	TEST_EQUAL(e.size(), 5)

	f.getOptions().setMZRange(makeRange(1025.0, 2000.0));
	f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"),e);
	TEST_EQUAL(e.size(), 3)

	f.getOptions().setIntensityRange(makeRange(290.0, 310.0));
	f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"),e);
	TEST_EQUAL(e.size(), 1)

	f.getOptions().setMetadataOnly(true);
	f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML"),e);
	TEST_EQUAL(e.getIdentifier(), "lsid2")
	TEST_EQUAL(e.size(), 0)
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
  FeatureXMLFile f;
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML")),true);
	TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_2_options.featureXML")),true);

	FeatureMap<> e;
	String filename;

  //test if empty file is valid
	NEW_TMP_FILE(filename)
	f.store(filename,e);
  TEST_EQUAL(f.isValid(filename),true);

	//test if full file is valid
	NEW_TMP_FILE(filename);
	f.load(OPENMS_GET_TEST_DATA_PATH("FeatureXMLFile_1.featureXML"),e);
	f.store(filename, e);
  TEST_EQUAL(f.isValid(filename),true);
END_SECTION

START_SECTION((const PeakFileOptions& getOptions() const))
	FeatureXMLFile f;
 	FeatureMap<> e;
	f.getOptions().setRTRange(makeRange(1.5, 4.5));
	f.getOptions().setIntensityRange(makeRange(290.0, 310.0));

	const PeakFileOptions pfo = f.getOptions();

	TEST_EQUAL(pfo.getRTRange(),makeRange(1.5, 4.5))
	TEST_EQUAL(pfo.getIntensityRange(),makeRange(290.0, 310.0))
END_SECTION

START_SECTION([EXTRA])
{
  Feature f1;
  f1.setRT(1001);
  f1.setMZ(1002);
  f1.setCharge(1003);
  Feature f1_cpy(f1);
  Feature f11;
  f11.setRT(1101);
  f11.setMZ(1102);
  Feature f12;
  f12.setRT(1201);
  f12.setMZ(1202);
  Feature f13;
  f13.setRT(1301);
  f13.setMZ(1302);
  TEST_EQUAL(f1.getSubordinates().empty(),true);
  f1.getSubordinates().push_back(f11);
  TEST_EQUAL(f1.getSubordinates().size(),1);
  f1.getSubordinates().push_back(f12);
  TEST_EQUAL(f1.getSubordinates().size(),2);
  f1.getSubordinates().push_back(f13);
  TEST_EQUAL(f1.getSubordinates().size(),3);
  TEST_EQUAL(f1.getRT(),1001);
  TEST_EQUAL(f1.getSubordinates()[0].getRT(),1101);
  TEST_EQUAL(f1.getSubordinates()[1].getRT(),1201);
  TEST_EQUAL(f1.getSubordinates()[2].getRT(),1301);
  const Feature &f1_cref = f1;
  TEST_EQUAL(f1_cref.getMZ(),1002);
  TEST_EQUAL(f1_cref.getSubordinates()[0].getMZ(),1102);
  TEST_EQUAL(f1_cref.getSubordinates()[1].getMZ(),1202);
  TEST_EQUAL(f1_cref.getSubordinates()[2].getMZ(),1302);
  TEST_NOT_EQUAL(f1_cref,f1_cpy);
  Feature f1_cpy2(f1);
  TEST_EQUAL(f1_cpy2,f1);
  f1.getSubordinates().clear();
  TEST_EQUAL(f1_cref,f1_cpy);

  Feature f2;
  f2.setRT(1001);
  f2.setMZ(1002);
  f2.setCharge(1003);
  TEST_NOT_EQUAL(f1_cpy2.getSubordinates().empty(),true);
  f2.setSubordinates(f1_cpy2.getSubordinates());
  TEST_EQUAL(f2,f1_cpy2);

  String filename;
  NEW_TMP_FILE(filename);
  FeatureXMLFile f;
  FeatureMap<> e;
  e.push_back(f1);
  e.push_back(f2);

  // this will print the number of newly assigned unique ids
  STATUS(e.applyMemberFunction(&UniqueIdInterface::ensureUniqueId));

  f.store(filename, e);
  FeatureMap<> e2;
  f.load(filename,e2);
  TEST_EQUAL(e==e2,true);
  String filename2;
  NEW_TMP_FILE(filename2);
  f.store(filename2, e2);

}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

