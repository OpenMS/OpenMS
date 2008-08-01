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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
///////////////////////////

#include <OpenMS/FORMAT/HANDLERS/MzDataExpSettHandler.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/Sample.h>

#include <sstream>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>


///////////////////////////

START_TEST(MzDataExpSettHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS; 
using namespace OpenMS::Internal; 
using namespace std; 
using namespace xercesc;

CHECK((MzDataExpSettHandler(ExperimentalSettings& exp, const String& filename)))
	XMLPlatformUtils::Initialize();
	SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
	parser->setFeature(XMLUni::fgSAX2CoreNameSpaces,false);
	parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes,false);
	
	ExperimentalSettings e;
	MzDataExpSettHandler handler(e,"");
	
	parser->setContentHandler(&handler);
	parser->setErrorHandler(&handler);
	
	LocalFileInputSource source(Internal::StringManager().convert("data/MzDataExpSett_test_1.xml") );
	parser->parse(source);
	delete(parser);

	PRECISION(0.01)

	//---------------------------------------------------------------------------
  // const SourceFile& getSourceFile() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSourceFile().getNameOfFile(), "MzDataFile_test_1.raw");
  TEST_EQUAL(e.getSourceFile().getPathToFile(), "/share/data/");
  TEST_EQUAL(e.getSourceFile().getFileType(), "MS");

  //---------------------------------------------------------------------------
  // const std::vector<ContactPerson>& getContacts() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getContacts().size(), 2);
  ABORT_IF(e.getContacts().size()!=2);
  TEST_EQUAL(e.getContacts()[0].getFirstName(), "John");
  TEST_EQUAL(e.getContacts()[0].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[0].getInstitution(), "department 1");
  TEST_EQUAL(e.getContacts()[0].getContactInfo(), "www.john.doe");
  TEST_EQUAL(e.getContacts()[1].getFirstName(), "Jane");
  TEST_EQUAL(e.getContacts()[1].getLastName(), "Doe");
  TEST_EQUAL(e.getContacts()[1].getInstitution(), "department 2");
  TEST_EQUAL(e.getContacts()[1].getContactInfo(), "www.jane.doe");

  //---------------------------------------------------------------------------
  // const Software& getSoftware() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getSoftware().getName(), "MS-X");
  TEST_EQUAL(e.getSoftware().getVersion(), "1.0");
  TEST_EQUAL(e.getSoftware().getComment(), "none");
	String tmp;
	e.getSoftware().getCompletionTime().get(tmp);
  TEST_EQUAL(tmp, "2001-02-03 04:05:06");

  //---------------------------------------------------------------------------
  // const ProcessingMethod& getProcessingMethod() const;
  //---------------------------------------------------------------------------
  TEST_EQUAL(e.getProcessingMethod().getDeisotoping(), false)
	TEST_EQUAL(e.getProcessingMethod().getChargeDeconvolution(), false)
	TEST_EQUAL(e.getProcessingMethod().getSpectrumType(), SpectrumSettings::PEAKS)
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e.getProcessingMethod().getMetaValue("ProcessingComment"), "Processed")

  //---------------------------------------------------------------------------
  // const Instrument& getInstrument() const;
  //---------------------------------------------------------------------------
	const Instrument& inst = e.getInstrument();
  TEST_EQUAL(inst.getName(), "MS-Instrument")
	TEST_EQUAL(inst.getVendor(), "MS-Vendor")
	TEST_EQUAL(inst.getModel(), "MS 1")
	TEST_EQUAL(inst.getCustomizations(), "tuned")
	TEST_EQUAL(inst.getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMetaValue("AdditionalComment"), "Additional")
	TEST_EQUAL(inst.getIonSource().getIonizationMethod(), IonSource::ESI)
	TEST_EQUAL(inst.getIonSource().getInletType(), IonSource::DIRECT)
	TEST_EQUAL(inst.getIonSource().getPolarity(), IonSource::NEGATIVE) 
	TEST_EQUAL(inst.getIonSource().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getIonSource().getMetaValue("SourceComment"), "Source")
	TEST_EQUAL(inst.getIonDetector().getType(), IonDetector::FARADAYCUP)
	TEST_EQUAL(inst.getIonDetector().getAcquisitionMode(), IonDetector::TDC)
	TEST_EQUAL(inst.getIonDetector().getResolution(), 0.815f)
	TEST_EQUAL(inst.getIonDetector().getADCSamplingFrequency(), 11.22f)
	TEST_EQUAL(inst.getIonDetector().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getIonDetector().getMetaValue("DetectorComment"), "Detector")
  TEST_EQUAL(inst.getMassAnalyzers().size(), 2)
  ABORT_IF(inst.getMassAnalyzers().size()!=2);
	TEST_EQUAL(inst.getMassAnalyzers()[0].getType(), MassAnalyzer::PAULIONTRAP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionMethod(), MassAnalyzer::FWHM)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolutionType(), MassAnalyzer::CONSTANT)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanFunction(), MassAnalyzer::MASSSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanDirection(), MassAnalyzer::UP)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanLaw(), MassAnalyzer::LINEAR)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTandemScanMethod(), MassAnalyzer::PRECURSORIONSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getReflectronState(), MassAnalyzer::OFF)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getResolution(), 22.33f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getAccuracy(), 33.44f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanRate(), 44.55f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getScanTime(), 55.66f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getTOFTotalPathLength(), 66.77f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getIsolationWidth(), 77.88f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getFinalMSExponent(), 2)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMagneticFieldStrength(), 88.99f)
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMassAnalyzers()[0].getMetaValue("AnalyzerComment"), "Analyzer 1")
	TEST_EQUAL(inst.getMassAnalyzers()[1].getType(), MassAnalyzer::QUADRUPOLE)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionMethod(), MassAnalyzer::BASELINE)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolutionType(), MassAnalyzer::PROPORTIONAL)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanFunction(), MassAnalyzer::SELECTEDIONDETECTION)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanDirection(), MassAnalyzer::DOWN)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanLaw(), MassAnalyzer::EXPONENTIAL)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getTandemScanMethod(), MassAnalyzer::PRODUCTIONSCAN)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getReflectronState(), MassAnalyzer::ON)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getResolution(), 12.3f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getAccuracy(), 13.4f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanRate(), 14.5f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getScanTime(), 15.6f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getTOFTotalPathLength(), 16.7f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getIsolationWidth(), 17.8f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getFinalMSExponent(), -2)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMagneticFieldStrength(), 18.9f)
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(inst.getMassAnalyzers()[1].getMetaValue("AnalyzerComment"), "Analyzer 2")

  //---------------------------------------------------------------------------
	// const Sample& getSample()
  //---------------------------------------------------------------------------
	TEST_EQUAL(e.getSample().getName(), "MS-Sample")
	TEST_EQUAL(e.getSample().getNumber(), "0-815")
	TEST_EQUAL(e.getSample().getState(), Sample::GAS)
 	TEST_EQUAL(e.getSample().getMass(), 1.01f)
	TEST_EQUAL(e.getSample().getVolume(), 2.02f)
	TEST_EQUAL(e.getSample().getConcentration(), 3.03f)
	TEST_EQUAL(e.getSample().getMetaValue("URL"), "www.open-ms.de")
	TEST_EQUAL(e.getSample().getMetaValue("SampleComment"), "Sample")
RESULT

CHECK((~MzDataExpSettHandler()))
	// implicitly tested in the test above
	NOT_TESTABLE
RESULT

CHECK((void characters(const XMLCh* const chars, unsigned int length)))
	// implicitly tested in the test above
	NOT_TESTABLE
RESULT

CHECK((void endElement(const XMLCh *const , const XMLCh *const, const XMLCh *const qname)))
	// implicitly tested in the test above
	NOT_TESTABLE
RESULT

CHECK((virtual void startElement(const XMLCh *const , const XMLCh *const, const XMLCh *const qname, const xercesc::Attributes &attributes)))
	// implicitly tested in the test above
	NOT_TESTABLE
RESULT

CHECK((MzDataExpSettHandler(const ExperimentalSettings& exp, const String& filename)))
	XMLPlatformUtils::Initialize();
	SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
	parser->setFeature(XMLUni::fgSAX2CoreNameSpaces,false);
	parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes,false);
	
	ExperimentalSettings e1;
	MzDataExpSettHandler handler(e1,"");
	
	parser->setContentHandler(&handler);
	parser->setErrorHandler(&handler);
	
	LocalFileInputSource source(Internal::StringManager().convert("data/MzDataExpSett_test_1.xml") );
	parser->parse(source);
	delete(parser);

	std::stringstream sstr;
	//read data and close stream
	const ExperimentalSettings& e2 = e1;
	Internal::MzDataExpSettHandler handler2(e2,"");
	handler2.writeTo(sstr);

	std::ifstream file("data/MzDataExpSett_test_1.xml");
	String s,t;
  while (getline(file,s))
  {
		getline(sstr,t);
		TEST_EQUAL(s.trim(),t.trim())
	}
RESULT

CHECK((void writeTo(std::ostream& os)))
	// implicitly tested in the test above
  NOT_TESTABLE	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
