// -*- Mode: C++; tab-width: 2; -*-
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

#include <iostream>

#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(PILISModel_test.C, "$Id: PILISModel_test.C 6446 2009-11-20 16:21:41Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISModel* ptr = 0;
const AASequence peptide("DFPIANGER");
START_SECTION(PILISModel())
	ptr = new PILISModel();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~PILISModel())
	delete ptr;
END_SECTION

ptr = new PILISModel();

START_SECTION(PILISModel(const PILISModel& model))
	PILISModel copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PILISModel& operator = (const PILISModel& mode))
	PILISModel copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void writeGraphMLFile(const String& filename))
	NOT_TESTABLE // will be tested in the next section, to avoid time consuming instantiation
END_SECTION

START_SECTION(void writeToFile(const String& filename))
	String filename;
	NEW_TMP_FILE(filename)
  PILISModel model;
  Param p(model.getParameters());
  p.setValue("model_depth", 2);
  p.setValue("visible_model_depth", 3);
  model.setParameters(p);
  model.init(true);
	model.writeToFile(filename);

	String graphml_filename;
	NEW_TMP_FILE(graphml_filename)
	model.writeGraphMLFile(graphml_filename);
	//TEST_FILE_SIMILAR(graphml_filename, OPENMS_GET_TEST_DATA_PATH("PILISModel_test.graphML"))

	PILISModel model2;
	model2.readFromFile(filename);
	TEST_EQUAL(model.getParameters() == model2.getParameters(), true)
	//model.getParameters().store("p.param");
	//model2.getParameters().store("p2.param");
	RichPeakSpectrum spec, spec2;
	model.getSpectrum(spec, "DER", 1);
	model2.getSpectrum(spec2, "DER", 1);
	TOLERANCE_ABSOLUTE(0.1) // arithmetic with small numbers...
	TEST_EQUAL(spec.size(), spec2.size())
	for (Size i = 0; i != spec.size(); ++i)
	{
		TEST_REAL_SIMILAR(spec[i].getMZ(), spec2[i].getMZ())
		TEST_REAL_SIMILAR(spec[i].getIntensity(), spec2[i].getIntensity())
	}
END_SECTION

START_SECTION(void readFromFile(const String& filename))
	NOT_TESTABLE // tested in previous section
END_SECTION

START_SECTION(void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, UInt charge))
	RichPeakSpectrum spec;
	PILISModel model;
	model.init(true);
	model.getSpectrum(spec, "DER", 1);
	TEST_EQUAL(spec.size(), 18)
END_SECTION

START_SECTION(void train(const RichPeakSpectrum&, const AASequence& peptide, UInt charge))
	RichPeakSpectrum spec1, spec2, spec3;
	PILISModel model;
	model.init(true);
	model.getSpectrum(spec1, "DER", 1);
	model.getSpectrum(spec2, "DEK", 1);
	model.getSpectrum(spec3, "DQK", 1);
	TEST_EQUAL(spec1.size(), 18);
	TEST_EQUAL(spec2.size(), 19);
	TEST_EQUAL(spec3.size(), 21);
	model.train(spec1, "DER", 1);
	model.train(spec2, "DEK", 1);
	model.train(spec3, "DQK", 1);

	// before calling evaluate the generated spectra should not change!
	RichPeakSpectrum spec4, spec5, spec6;
	model.getSpectrum(spec4, "DER", 1);
	model.getSpectrum(spec5, "DEK", 1);
	model.getSpectrum(spec6, "DQK", 1);
	TEST_EQUAL(spec1.size(), spec4.size())
	TEST_EQUAL(spec2.size(), spec5.size())
	TEST_EQUAL(spec3.size(), spec6.size())
  for (Size i = 0; i != spec1.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec1[i].getMZ(), spec4[i].getMZ())
    TEST_REAL_SIMILAR(spec1[i].getIntensity(), spec4[i].getIntensity())
  }
  for (Size i = 0; i != spec2.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec2[i].getMZ(), spec5[i].getMZ())
    TEST_REAL_SIMILAR(spec2[i].getIntensity(), spec5[i].getIntensity())
  }
  for (Size i = 0; i != spec3.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec3[i].getMZ(), spec6[i].getMZ())
    TEST_REAL_SIMILAR(spec3[i].getIntensity(), spec6[i].getIntensity())
  }

	model.evaluate();

	RichPeakSpectrum spec;
	model.getSpectrum(spec, "DER", 1);

/*	
	for (Size i = 0; i != spec.size(); ++i)
	{
		cout << i << " " << spec[i].getMZ() << " " << spec[i].getIntensity() << " " << spec[i].getMetaValue("IonName") << endl;
	}

	for (Size i = 0; i != spec1.size(); ++i)
	{
		cout << i << " " << spec1[i].getMZ() << " " << spec1[i].getIntensity() << " " << spec1[i].getMetaValue("IonName") << endl;
	}
*/

	TEST_NOT_EQUAL(spec == spec1, true)
	

END_SECTION

START_SECTION(void evaluate())
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void init(bool generate_models=true)))
	NOT_TESTABLE // tested implicetely above
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

