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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/MassAnalyzer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassAnalyzer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassAnalyzer* ptr = 0;
CHECK(MassAnalyzer())
	ptr = new MassAnalyzer();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MassAnalyzer())
	delete ptr;
RESULT

CHECK(AnalyzerType getType() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getType(),MassAnalyzer::ANALYZERNULL);
RESULT

CHECK(ReflectronState getReflectronState() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::REFLSTATENULL);
RESULT

CHECK(ResolutionMethod getResolutionMethod() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::RESMETHNULL);
RESULT

CHECK(ResolutionType getResolutionType() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::RESTYPENULL);
RESULT

CHECK(ScanDirection getScanDirection() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::SCANDIRNULL);
RESULT

CHECK(ScanFunction getScanFunction() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanFunction(),MassAnalyzer::SCANFCTNULL);
RESULT

CHECK(ScanLaw getScanLaw() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::SCANLAWNULL);
RESULT

CHECK(SignedInt getFinalMSExponent() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getFinalMSExponent(),0);
RESULT

CHECK(TandemScanningMethod getTandemScanMethod() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getTandemScanMethod(),MassAnalyzer::TANDEMNULL);
RESULT

CHECK(float getAccuracy() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getAccuracy(),0.0);
RESULT

CHECK(float getIsolationWidth() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getIsolationWidth(),0.0);
RESULT

CHECK(float getMagneticFieldStrength() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getMagneticFieldStrength(),0.0);
RESULT

CHECK(float getResolution() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getResolution(),0.0);
RESULT

CHECK(float getScanRate() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getScanRate(),0.0);
RESULT

CHECK(float getScanTime() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getScanTime(),0.0);
RESULT

CHECK(float getTOFTotalPathLength() const)
  MassAnalyzer tmp;
  TEST_REAL_EQUAL(tmp.getTOFTotalPathLength(),0.0);
RESULT

CHECK(void setType(AnalyzerType type))
  MassAnalyzer tmp;
  tmp.setType(MassAnalyzer::QUADRUPOLE);
  TEST_EQUAL(tmp.getType(),MassAnalyzer::QUADRUPOLE);
RESULT

CHECK(void setAccuracy(float accuracy))
  MassAnalyzer tmp;
  tmp.setAccuracy(47.11);
  TEST_REAL_EQUAL(tmp.getAccuracy(),47.11);
RESULT

CHECK(void setFinalMSExponent(SignedInt final_MS_exponent))
  MassAnalyzer tmp;
  tmp.setFinalMSExponent(47);
  TEST_EQUAL(tmp.getFinalMSExponent(),47);
RESULT

CHECK(void setIsolationWidth(float isolation_width))
  MassAnalyzer tmp;
  tmp.setIsolationWidth(47.11);
  TEST_REAL_EQUAL(tmp.getIsolationWidth(),47.11);
RESULT

CHECK(void setMagneticFieldStrength(float magnetic_field_strength))
  MassAnalyzer tmp;
  tmp.setMagneticFieldStrength(47.11);
  TEST_REAL_EQUAL(tmp.getMagneticFieldStrength(),47.11);
RESULT

CHECK(void setReflectronState(ReflectronState reflecton_state))
  MassAnalyzer tmp;
  tmp.setReflectronState(MassAnalyzer::ON);
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ON);
RESULT

CHECK(void setResolution(float resolution))
  MassAnalyzer tmp;
  tmp.setResolution(47.11);
  TEST_REAL_EQUAL(tmp.getResolution(),47.11);
RESULT

CHECK(void setResolutionMethod(ResolutionMethod resolution_method))
  MassAnalyzer tmp;
  tmp.setResolutionMethod(MassAnalyzer::FWHM);
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::FWHM);
RESULT

CHECK(void setResolutionType(ResolutionType resolution_type))
  MassAnalyzer tmp;
  tmp.setResolutionType(MassAnalyzer::CONSTANT);
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::CONSTANT);
RESULT

CHECK(void setScanDirection(ScanDirection scan_direction))
  MassAnalyzer tmp;
  tmp.setScanDirection(MassAnalyzer::UP);
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::UP);
RESULT

CHECK(void setScanFunction(ScanFunction scan_function))
  MassAnalyzer tmp;
  tmp.setScanFunction(MassAnalyzer::MASSSCAN);
  TEST_EQUAL(tmp.getScanFunction(),MassAnalyzer::MASSSCAN);
RESULT

CHECK(void setScanLaw(ScanLaw scan_law))
  MassAnalyzer tmp;
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::LINEAR);
RESULT

CHECK(void setScanRate(float scan_rate))
  MassAnalyzer tmp;
  tmp.setScanRate(47.11);
  TEST_REAL_EQUAL(tmp.getScanRate(),47.11);
RESULT

CHECK(void setScanTime(float scan_time))
  MassAnalyzer tmp;
  tmp.setScanTime(47.11);
  TEST_REAL_EQUAL(tmp.getScanTime(),47.11);
RESULT

CHECK(void setTOFTotalPathLength(float TOF_total_path_length))
  MassAnalyzer tmp;
  tmp.setTOFTotalPathLength(47.11);
  TEST_REAL_EQUAL(tmp.getTOFTotalPathLength(),47.11);
RESULT

CHECK(void setTandemScanMethod(TandemScanningMethod tandem_scan_method))
  MassAnalyzer tmp;
  tmp.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
  TEST_EQUAL(tmp.getTandemScanMethod(),MassAnalyzer::PRODUCTIONSCAN);
RESULT

CHECK(MassAnalyzer(const MassAnalyzer& source))
  MassAnalyzer tmp;
  tmp.setType(MassAnalyzer::QUADRUPOLE);
  tmp.setAccuracy(47.11);
  tmp.setFinalMSExponent(47);
  tmp.setIsolationWidth(47.12);
  tmp.setMagneticFieldStrength(47.13);
  tmp.setReflectronState(MassAnalyzer::ON);
  tmp.setResolution(47.14);
  tmp.setResolutionMethod(MassAnalyzer::FWHM);
  tmp.setResolutionType(MassAnalyzer::CONSTANT);
  tmp.setScanDirection(MassAnalyzer::UP);
  tmp.setScanFunction(MassAnalyzer::MASSSCAN);
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  tmp.setScanRate(47.15);
  tmp.setScanTime(47.16);
  tmp.setTOFTotalPathLength(47.17);
  tmp.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
	tmp.setMetaValue("label",String("label"));
	
	MassAnalyzer tmp2(tmp);
	TEST_EQUAL(tmp.getType(),MassAnalyzer::QUADRUPOLE);
	TEST_REAL_EQUAL(tmp.getAccuracy(),47.11);
	TEST_EQUAL(tmp.getFinalMSExponent(),47);
	TEST_REAL_EQUAL(tmp.getIsolationWidth(),47.12);
	TEST_REAL_EQUAL(tmp.getMagneticFieldStrength(),47.13);
	TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ON);
	TEST_REAL_EQUAL(tmp.getResolution(),47.14);
	TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::FWHM);
	TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::CONSTANT);
	TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::UP);
	TEST_EQUAL(tmp.getScanFunction(),MassAnalyzer::MASSSCAN);
	TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::LINEAR);
	TEST_REAL_EQUAL(tmp.getScanRate(),47.15);
	TEST_REAL_EQUAL(tmp.getScanTime(),47.16);
	TEST_REAL_EQUAL(tmp.getTOFTotalPathLength(),47.17);
	TEST_EQUAL(tmp.getTandemScanMethod(),MassAnalyzer::PRODUCTIONSCAN);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
RESULT

CHECK(MassAnalyzer& operator= (const MassAnalyzer& source))
  MassAnalyzer tmp;
  
  tmp.setType(MassAnalyzer::QUADRUPOLE);
  tmp.setAccuracy(47.11);
  tmp.setFinalMSExponent(47);
  tmp.setIsolationWidth(47.12);
  tmp.setMagneticFieldStrength(47.13);
  tmp.setReflectronState(MassAnalyzer::ON);
  tmp.setResolution(47.14);
  tmp.setResolutionMethod(MassAnalyzer::FWHM);
  tmp.setResolutionType(MassAnalyzer::CONSTANT);
  tmp.setScanDirection(MassAnalyzer::UP);
  tmp.setScanFunction(MassAnalyzer::MASSSCAN);
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  tmp.setScanRate(47.15);
  tmp.setScanTime(47.16);
  tmp.setTOFTotalPathLength(47.17);
  tmp.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
	tmp.setMetaValue("label",String("label"));
	
	MassAnalyzer tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getType(),MassAnalyzer::QUADRUPOLE);
	TEST_REAL_EQUAL(tmp2.getAccuracy(),47.11);
	TEST_EQUAL(tmp2.getFinalMSExponent(),47);
	TEST_REAL_EQUAL(tmp2.getIsolationWidth(),47.12);
	TEST_REAL_EQUAL(tmp2.getMagneticFieldStrength(),47.13);
	TEST_EQUAL(tmp2.getReflectronState(),MassAnalyzer::ON);
	TEST_REAL_EQUAL(tmp2.getResolution(),47.14);
	TEST_EQUAL(tmp2.getResolutionMethod(),MassAnalyzer::FWHM);
	TEST_EQUAL(tmp2.getResolutionType(),MassAnalyzer::CONSTANT);
	TEST_EQUAL(tmp2.getScanDirection(),MassAnalyzer::UP);
	TEST_EQUAL(tmp2.getScanFunction(),MassAnalyzer::MASSSCAN);
	TEST_EQUAL(tmp2.getScanLaw(),MassAnalyzer::LINEAR);
	TEST_REAL_EQUAL(tmp2.getScanRate(),47.15);
	TEST_REAL_EQUAL(tmp2.getScanTime(),47.16);
	TEST_REAL_EQUAL(tmp2.getTOFTotalPathLength(),47.17);
	TEST_EQUAL(tmp2.getTandemScanMethod(),MassAnalyzer::PRODUCTIONSCAN);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");

	tmp2 = MassAnalyzer();
	TEST_EQUAL(tmp2.getType(),MassAnalyzer::ANALYZERNULL);
	TEST_REAL_EQUAL(tmp2.getAccuracy(),0.0);
	TEST_EQUAL(tmp2.getFinalMSExponent(),0);
	TEST_REAL_EQUAL(tmp2.getIsolationWidth(),0.0);
	TEST_REAL_EQUAL(tmp2.getMagneticFieldStrength(),0.0);
	TEST_EQUAL(tmp2.getReflectronState(),MassAnalyzer::REFLSTATENULL);
	TEST_REAL_EQUAL(tmp2.getResolution(),0.0);
	TEST_EQUAL(tmp2.getResolutionMethod(),MassAnalyzer::RESMETHNULL);
	TEST_EQUAL(tmp2.getResolutionType(),MassAnalyzer::RESTYPENULL);
	TEST_EQUAL(tmp2.getScanDirection(),MassAnalyzer::SCANDIRNULL);
	TEST_EQUAL(tmp2.getScanFunction(),MassAnalyzer::SCANFCTNULL);
	TEST_EQUAL(tmp2.getScanLaw(),MassAnalyzer::SCANLAWNULL);
	TEST_REAL_EQUAL(tmp2.getScanRate(),0.0);
	TEST_REAL_EQUAL(tmp2.getScanTime(),0.0);
	TEST_REAL_EQUAL(tmp2.getTOFTotalPathLength(),0.0);
	TEST_EQUAL(tmp2.getTandemScanMethod(),MassAnalyzer::TANDEMNULL);
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK(bool operator== (const MassAnalyzer& rhs) const)
  MassAnalyzer edit, empty;
  
	TEST_EQUAL(edit==empty,true);

	edit=empty;
	edit.setType(MassAnalyzer::QUADRUPOLE);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setAccuracy(47.11);
	TEST_EQUAL(edit==empty,false);

	edit=empty;
	edit.setFinalMSExponent(47);
	TEST_EQUAL(edit==empty,false);

	edit=empty;
	edit.setIsolationWidth(47.12);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setMagneticFieldStrength(47.13);
	TEST_EQUAL(edit==empty,false);	
	
	edit=empty;
	edit.setReflectronState(MassAnalyzer::ON);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setResolution(47.14);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setResolutionMethod(MassAnalyzer::FWHM);
	TEST_EQUAL(edit==empty,false);

	edit=empty;
	edit.setResolutionType(MassAnalyzer::CONSTANT);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanDirection(MassAnalyzer::UP);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanFunction(MassAnalyzer::MASSSCAN);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanLaw(MassAnalyzer::LINEAR);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanRate(47.15);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanTime(47.16);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setTOFTotalPathLength(47.17);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
	edit.setMetaValue("label",String("label"));
RESULT

CHECK(bool operator!= (const MassAnalyzer& rhs) const)
  MassAnalyzer edit, empty;
  
	TEST_EQUAL(edit!=empty,false);

	edit=empty;
	edit.setType(MassAnalyzer::QUADRUPOLE);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setAccuracy(47.11);
	TEST_EQUAL(edit!=empty,true);

	edit=empty;
	edit.setFinalMSExponent(47);
	TEST_EQUAL(edit!=empty,true);

	edit=empty;
	edit.setIsolationWidth(47.12);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setMagneticFieldStrength(47.13);
	TEST_EQUAL(edit!=empty,true);	
	
	edit=empty;
	edit.setReflectronState(MassAnalyzer::ON);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setResolution(47.14);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setResolutionMethod(MassAnalyzer::FWHM);
	TEST_EQUAL(edit!=empty,true);

	edit=empty;
	edit.setResolutionType(MassAnalyzer::CONSTANT);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanDirection(MassAnalyzer::UP);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanFunction(MassAnalyzer::MASSSCAN);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanLaw(MassAnalyzer::LINEAR);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanRate(47.15);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanTime(47.16);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setTOFTotalPathLength(47.17);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
	edit.setMetaValue("label",String("label"));
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



