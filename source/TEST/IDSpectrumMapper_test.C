// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/ID/IDSpectrumMapper.h>
#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DSpectrum.h>

///////////////////////////

START_TEST(IDSpectrumMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IDSpectrumMapper annotator;
IDSpectrumMapper* ptr;
vector<IdentificationData> identifications; 
vector<ProteinIdentification> protein_identifications;
float precision = 0.1;
CHECK((IDSpectrumMapper()))
	ptr = new IDSpectrumMapper();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((template<class PeakT> UnsignedInt annotate(MSExperiment< PeakT >& experiment, const std::vector<IdentificationData>& identifications, float precision = 0.01f)))
	vector<IdentificationData> identifications2; 
	MSExperiment< DPeak<1> > experiment;
	MSSpectrum< DPeak<1> > spectrum;
	DSpectrum< 1 >::PrecursorPeakType peak;
	AnalysisXMLFile().load("data/IDSpectrumMapper_test.analysisXML",
								protein_identifications, 
					   		identifications);
	
	peak = spectrum.getPrecursorPeak();
	peak.setPosition(0);
	spectrum.setRetentionTime(60);
	experiment.push_back(spectrum);							
	experiment[0].setPrecursorPeak(peak);
	peak.setPosition(20);
	spectrum.setRetentionTime(181);
	experiment.push_back(spectrum);							
	experiment[1].setPrecursorPeak(peak);
	peak.setPosition(11);
	spectrum.setRetentionTime(120.0001);
	experiment.push_back(spectrum);							
	experiment[2].setPrecursorPeak(peak);
	

	annotator.annotate(experiment, identifications, precision);
  annotator.getAnnotations(experiment, identifications2);
  TEST_EQUAL(identifications2.size(), 2)
  PRECISION(precision);
  TEST_REAL_EQUAL(identifications2[0].rt, identifications[0].rt)
  TEST_REAL_EQUAL(identifications2[0].mz, identifications[0].mz)
  TEST_REAL_EQUAL(identifications2[1].rt, identifications[1].rt)
  TEST_REAL_EQUAL(identifications2[1].mz, identifications[1].mz)
RESULT

CHECK((template<class PeakT> void getAnnotations(const MSExperiment< PeakT >& experiment, std::vector<IdentificationData>& identifications)))
	vector<IdentificationData> identifications2; 
	MSExperiment< DPeak<1> > experiment;
	MSSpectrum< DPeak<1> > spectrum;
	DSpectrum< 1 >::PrecursorPeakType peak;
	AnalysisXMLFile().load("data/IDSpectrumMapper_test.analysisXML",
								protein_identifications, 
					   		identifications);
	
	peak = spectrum.getPrecursorPeak();
	peak.setPosition(0);
	spectrum.setRetentionTime(60);
	experiment.push_back(spectrum);							
	experiment[0].setPrecursorPeak(peak);
	peak.setPosition(11);
	spectrum.setRetentionTime(120);
	experiment.push_back(spectrum);							
	experiment[1].setPrecursorPeak(peak);
	peak.setPosition(20);
	spectrum.setRetentionTime(180);
	experiment.push_back(spectrum);							
	experiment[2].setPrecursorPeak(peak);
	
	annotator.annotate(experiment, identifications, precision);
  annotator.getAnnotations(experiment, identifications2);
  TEST_EQUAL(identifications2.size(), identifications.size())
  PRECISION(precision);
  TEST_REAL_EQUAL(identifications2[0].rt, identifications[0].rt)
  TEST_REAL_EQUAL(identifications2[0].mz, identifications[0].mz)
  TEST_REAL_EQUAL(identifications2[1].rt, identifications[1].rt)
  TEST_REAL_EQUAL(identifications2[1].mz, identifications[1].mz)
  TEST_REAL_EQUAL(identifications2[2].rt, identifications[2].rt)
  TEST_REAL_EQUAL(identifications2[2].mz, identifications[2].mz)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
