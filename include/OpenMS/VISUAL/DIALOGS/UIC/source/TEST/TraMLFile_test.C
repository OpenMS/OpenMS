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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace OpenMS;
using namespace std;

START_TEST(TraMLFile, "$Id: TraMLFile_test.C 6504 2010-01-04 13:53:50Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


TraMLFile* ptr = 0;
START_SECTION((TraMLFile()))
	ptr = new TraMLFile;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~TraMLFile()))
	delete ptr;
END_SECTION

START_SECTION((void load(const String& filename, MRMexperiment& exp)))
	TraMLFile file;
	TargetedExperiment exp;
	file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp);

END_SECTION

START_SECTION((void store(const String& filename, const TargetedExperiment& exp) const))
	TraMLFile file;
	
	{
		//load map
		TargetedExperiment exp_original;
		file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"), exp_original);
	 	
		//store map
		std::string tmp_filename;
	 	NEW_TMP_FILE(tmp_filename);
		file.store(tmp_filename, exp_original);
		//file.store("include.traML", TargetedExperiment());
		//load written map
		TargetedExperiment exp;
		file.load(tmp_filename, exp);

		//test if everything worked
		//TEST_EQUAL(exp == exp_original,true)
	}
/*
	//test with empty map
	{
		
		MSExperiment<> empty, exp;
		
		//this will be set when writing (forced by mzML)
		empty.getInstrument().getSoftware().setName("custom unreleased software tool");
		
		std::string tmp_filename;
		NEW_TMP_FILE(tmp_filename);
		file.store(tmp_filename,empty);
		file.load(tmp_filename,exp);
		TEST_EQUAL(exp==empty,true)
	}
	
	//test with one empty spectrum
	{
		MSExperiment<> empty, exp;
		empty.resize(1);
		
		//this will be set when writing (forced by mzML)
		empty.getInstrument().getSoftware().setName("custom unreleased software tool");
		empty[0].setNativeID("spectrum=0");
		empty[0].getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);
		empty[0].getDataProcessing().resize(1);
		empty[0].getDataProcessing()[0].getProcessingActions().insert(DataProcessing::CONVERSION_MZML);
		empty[0].getDataProcessing()[0].getSoftware().setName("custom unreleased software tool");
		
		std::string tmp_filename;
		NEW_TMP_FILE(tmp_filename);
		file.store(tmp_filename,empty);
		file.load(tmp_filename,exp);
		TEST_EQUAL(exp==empty,true)

		TEST_EQUAL(exp.size()==empty.size(),true)
		TEST_EQUAL(exp.ExperimentalSettings::operator==(empty),true)
		TEST_EQUAL(exp[0].SpectrumSettings::operator==(empty[0]),true)
		TEST_EQUAL(exp[0]==empty[0],true);
	}
*/
END_SECTION

START_SECTION(bool isValid(const String& filename, std::ostream& os = std::cerr))
	std::string tmp_filename;
  TraMLFile file;
  TargetedExperiment e;

  //written empty file
	NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isValid(tmp_filename), true);

	//written filled file
	NEW_TMP_FILE(tmp_filename);
	file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"),e);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isValid(tmp_filename),true);

END_SECTION

START_SECTION(bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings))

	std::string tmp_filename;
	TraMLFile file;
	StringList errors, warnings;
  TargetedExperiment e;

  //written empty file
	NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename,e);
  TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)

	//written filled file
	NEW_TMP_FILE(tmp_filename);
	file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.TraML"),e);
  file.store(tmp_filename,e);
  //TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
	TEST_EQUAL(errors.size(),0)
	TEST_EQUAL(warnings.size(),0)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

