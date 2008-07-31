// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/DSpectrum.h>

///////////////////////////

START_TEST(IDSpectrumMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

IDSpectrumMapper* ptr;
CHECK((IDSpectrumMapper()))
	ptr = new IDSpectrumMapper();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(template <class PeakT> UInt annotate(MSExperiment< PeakT > &experiment, const std::vector< PeptideIdentification > &identifications, DoubleReal precision=0.01f) throw (Exception::MissingInformation))
	//load id
	vector<PeptideIdentification> identifications; 
	vector<ProteinIdentification> protein_identifications;
	IdXMLFile().load("data/IDSpectrumMapper_test.idXML",
								protein_identifications, 
					   		identifications);

	TEST_EQUAL(identifications[0].getHits().size(), 2)
	TEST_EQUAL(identifications[1].getHits().size(), 1)
	TEST_EQUAL(identifications[2].getHits().size(), 2)
	TEST_EQUAL(protein_identifications[0].getHits().size(), 2)

	//create experiment
	MSExperiment<> experiment;
	MSSpectrum<> spectrum;
	DSpectrum<>::PrecursorPeakType peak;	
	peak = spectrum.getPrecursorPeak();
	peak.setPosition(0);
	spectrum.setRT(60);
	experiment.push_back(spectrum);							
	experiment[0].setPrecursorPeak(peak);
	peak.setPosition(20);
	spectrum.setRT(181);
	experiment.push_back(spectrum);							
	experiment[1].setPrecursorPeak(peak);
	peak.setPosition(11);
	spectrum.setRT(120.0001);
	experiment.push_back(spectrum);							
	experiment[2].setPrecursorPeak(peak);
	
	//map
	IDSpectrumMapper().annotate(experiment, identifications, 0.1);

	//scan 1
	TEST_EQUAL(experiment[0].getPeptideIdentifications().size(), 1)
	TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits().size(), 2)
	TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), "LHASGITVTEIPVTATNFK")
	TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits()[1].getSequence(), "MRSLGYVAVISAVATDTDK")
	
	//scan 2
	TEST_EQUAL(experiment[1].getPeptideIdentifications().size(), 0)
	
	//scan 3
	TEST_EQUAL(experiment[2].getPeptideIdentifications().size(), 1)	
	TEST_EQUAL(experiment[2].getPeptideIdentifications()[0].getHits().size(), 1)
	TEST_EQUAL(experiment[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), "HSKLSAK")

	
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
