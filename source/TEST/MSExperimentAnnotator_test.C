// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/MSExperimentAnnotator.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DPosition.h>
#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(MSExperimentAnnotator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MSExperimentAnnotator annotator;
MSExperimentAnnotator* ptr;
MzXMLFile file;
MSExperiment< DPeak<1> > experiment;
vector<Identification> identifications; 
vector<ProteinIdentification> protein_identifications; 
AnalysisXMLFile xml_file;
vector<float> precursor_retention_times;
vector<float> precursor_mz_values;
ContactPerson contact_person;
float precision = 0.1;

file.load("data/MzXMLFile_test_3.mzXML", experiment);

xml_file.load("data/MSExperimentAnnotatorFile_test.analysisXML",
							&protein_identifications, 
				   		&identifications, 
							&precursor_retention_times, 
							&precursor_mz_values, 
							&contact_person);								

CHECK(MSExperimentAnnotator())
	ptr = new MSExperimentAnnotator();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((template<class PeakT> UnsignedInt annotate(MSExperiment< PeakT >& experiment, const std::vector<Identification>& identifications, const std::vector<float>& precursor_retention_times, const std::vector<float>& precursor_mz_values, float precision = 0.01f)))
	vector<Identification> identifications2; 
	vector<float> precursor_retention_times2;
	vector<float> precursor_mz_values2;

	annotator.annotate(experiment,
      					 	   identifications,
      				 			 precursor_retention_times,
      				 			 precursor_mz_values,
      				 			 precision);
  annotator.getAnnotations(experiment,
      					 	   			&identifications2,
      				 			 			&precursor_retention_times2,
      				 			 			&precursor_mz_values2);
  TEST_EQUAL(identifications.size(), identifications2.size())
  PRECISION(precision);
  TEST_REAL_EQUAL(precursor_retention_times[0], precursor_retention_times2[0])
  TEST_REAL_EQUAL(precursor_mz_values[0], precursor_mz_values2[0])
RESULT

CHECK((template<class PeakT> void getAnnotations(const MSExperiment< PeakT >& experiment, std::vector<Identification>* identifications, std::vector<float>* precursor_retention_times, std::vector<float>* precursor_mz_values)))
	vector<Identification> identifications3; 
	vector<float> precursor_retention_times3;
	vector<float> precursor_mz_values3;

  annotator.getAnnotations(experiment,
      					 	   			&identifications3,
      				 			 			&precursor_retention_times3,
      				 			 			&precursor_mz_values3);
  TEST_EQUAL(identifications.size(), identifications3.size())
  PRECISION(precision);
  TEST_REAL_EQUAL(precursor_retention_times[0], precursor_retention_times3[0])
  TEST_REAL_EQUAL(precursor_mz_values[0], precursor_mz_values3[0])
  TEST_REAL_EQUAL(60, precursor_retention_times3[0])
  TEST_REAL_EQUAL(precursor_mz_values[0], 0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
