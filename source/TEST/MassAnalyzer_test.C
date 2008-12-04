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
#include <OpenMS/METADATA/MassAnalyzer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassAnalyzer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassAnalyzer* ptr = 0;
START_SECTION(MassAnalyzer())
	ptr = new MassAnalyzer();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~MassAnalyzer())
	delete ptr;
END_SECTION

START_SECTION(Int getOrder() const)
	MassAnalyzer tmp;
	TEST_EQUAL(tmp.getOrder(),0)
END_SECTION

START_SECTION(void setOrder(Int order))
	MassAnalyzer tmp;
	tmp.setOrder(4711);
	TEST_EQUAL(tmp.getOrder(),4711)
END_SECTION

START_SECTION(AnalyzerType getType() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getType(),MassAnalyzer::ANALYZERNULL);
END_SECTION

START_SECTION(ReflectronState getReflectronState() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::REFLSTATENULL);
END_SECTION

START_SECTION(ResolutionMethod getResolutionMethod() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::RESMETHNULL);
END_SECTION

START_SECTION(ResolutionType getResolutionType() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::RESTYPENULL);
END_SECTION

START_SECTION(ScanDirection getScanDirection() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::SCANDIRNULL);
END_SECTION

START_SECTION(ScanLaw getScanLaw() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::SCANLAWNULL);
END_SECTION

START_SECTION(Int getFinalMSExponent() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getFinalMSExponent(),0);
END_SECTION

START_SECTION(TandemScanningMethod getTandemScanMethod() const)
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getTandemScanMethod(),MassAnalyzer::TANDEMNULL);
END_SECTION

START_SECTION(float getAccuracy() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getAccuracy(),0.0);
END_SECTION

START_SECTION(float getIsolationWidth() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWidth(),0.0);
END_SECTION

START_SECTION(float getMagneticFieldStrength() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getMagneticFieldStrength(),0.0);
END_SECTION

START_SECTION(float getResolution() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getResolution(),0.0);
END_SECTION

START_SECTION(float getScanRate() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getScanRate(),0.0);
END_SECTION

START_SECTION(float getScanTime() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getScanTime(),0.0);
END_SECTION

START_SECTION(float getTOFTotalPathLength() const)
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getTOFTotalPathLength(),0.0);
END_SECTION

START_SECTION(void setType(AnalyzerType type))
  MassAnalyzer tmp;
  tmp.setType(MassAnalyzer::QUADRUPOLE);
  TEST_EQUAL(tmp.getType(),MassAnalyzer::QUADRUPOLE);
END_SECTION

START_SECTION(void setAccuracy(float accuracy))
  MassAnalyzer tmp;
  tmp.setAccuracy(47.11);
  TEST_REAL_SIMILAR(tmp.getAccuracy(),47.11);
END_SECTION

START_SECTION(void setFinalMSExponent(Int final_MS_exponent))
  MassAnalyzer tmp;
  tmp.setFinalMSExponent(47);
  TEST_EQUAL(tmp.getFinalMSExponent(),47);
END_SECTION

START_SECTION(void setIsolationWidth(float isolation_width))
  MassAnalyzer tmp;
  tmp.setIsolationWidth(47.11);
  TEST_REAL_SIMILAR(tmp.getIsolationWidth(),47.11);
END_SECTION

START_SECTION(void setMagneticFieldStrength(float magnetic_field_strength))
  MassAnalyzer tmp;
  tmp.setMagneticFieldStrength(47.11);
  TEST_REAL_SIMILAR(tmp.getMagneticFieldStrength(),47.11);
END_SECTION

START_SECTION(void setReflectronState(ReflectronState reflecton_state))
  MassAnalyzer tmp;
  tmp.setReflectronState(MassAnalyzer::ON);
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ON);
END_SECTION

START_SECTION(void setResolution(float resolution))
  MassAnalyzer tmp;
  tmp.setResolution(47.11);
  TEST_REAL_SIMILAR(tmp.getResolution(),47.11);
END_SECTION

START_SECTION(void setResolutionMethod(ResolutionMethod resolution_method))
  MassAnalyzer tmp;
  tmp.setResolutionMethod(MassAnalyzer::FWHM);
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::FWHM);
END_SECTION

START_SECTION(void setResolutionType(ResolutionType resolution_type))
  MassAnalyzer tmp;
  tmp.setResolutionType(MassAnalyzer::CONSTANT);
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::CONSTANT);
END_SECTION

START_SECTION(void setScanDirection(ScanDirection scan_direction))
  MassAnalyzer tmp;
  tmp.setScanDirection(MassAnalyzer::UP);
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::UP);
END_SECTION

START_SECTION(void setScanLaw(ScanLaw scan_law))
  MassAnalyzer tmp;
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::LINEAR);
END_SECTION

START_SECTION(void setScanRate(float scan_rate))
  MassAnalyzer tmp;
  tmp.setScanRate(47.11);
  TEST_REAL_SIMILAR(tmp.getScanRate(),47.11);
END_SECTION

START_SECTION(void setScanTime(float scan_time))
  MassAnalyzer tmp;
  tmp.setScanTime(47.11);
  TEST_REAL_SIMILAR(tmp.getScanTime(),47.11);
END_SECTION

START_SECTION(void setTOFTotalPathLength(float TOF_total_path_length))
  MassAnalyzer tmp;
  tmp.setTOFTotalPathLength(47.11);
  TEST_REAL_SIMILAR(tmp.getTOFTotalPathLength(),47.11);
END_SECTION

START_SECTION(void setTandemScanMethod(TandemScanningMethod tandem_scan_method))
  MassAnalyzer tmp;
  tmp.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
  TEST_EQUAL(tmp.getTandemScanMethod(),MassAnalyzer::PRODUCTIONSCAN);
END_SECTION

START_SECTION(MassAnalyzer(const MassAnalyzer& source))
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
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  tmp.setScanRate(47.15);
  tmp.setScanTime(47.16);
  tmp.setTOFTotalPathLength(47.17);
  tmp.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
	tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
	
	MassAnalyzer tmp2(tmp);
	TEST_EQUAL(tmp.getType(),MassAnalyzer::QUADRUPOLE);
	TEST_REAL_SIMILAR(tmp.getAccuracy(),47.11);
	TEST_EQUAL(tmp.getFinalMSExponent(),47);
	TEST_REAL_SIMILAR(tmp.getIsolationWidth(),47.12);
	TEST_REAL_SIMILAR(tmp.getMagneticFieldStrength(),47.13);
	TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ON);
	TEST_REAL_SIMILAR(tmp.getResolution(),47.14);
	TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::FWHM);
	TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::CONSTANT);
	TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::UP);
	TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::LINEAR);
	TEST_REAL_SIMILAR(tmp.getScanRate(),47.15);
	TEST_REAL_SIMILAR(tmp.getScanTime(),47.16);
	TEST_REAL_SIMILAR(tmp.getTOFTotalPathLength(),47.17);
	TEST_EQUAL(tmp.getTandemScanMethod(),MassAnalyzer::PRODUCTIONSCAN);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION(MassAnalyzer& operator= (const MassAnalyzer& source))
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
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  tmp.setScanRate(47.15);
  tmp.setScanTime(47.16);
  tmp.setTOFTotalPathLength(47.17);
  tmp.setTandemScanMethod(MassAnalyzer::PRODUCTIONSCAN);
	tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
	
	MassAnalyzer tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getType(),MassAnalyzer::QUADRUPOLE);
	TEST_REAL_SIMILAR(tmp2.getAccuracy(),47.11);
	TEST_EQUAL(tmp2.getFinalMSExponent(),47);
	TEST_REAL_SIMILAR(tmp2.getIsolationWidth(),47.12);
	TEST_REAL_SIMILAR(tmp2.getMagneticFieldStrength(),47.13);
	TEST_EQUAL(tmp2.getReflectronState(),MassAnalyzer::ON);
	TEST_REAL_SIMILAR(tmp2.getResolution(),47.14);
	TEST_EQUAL(tmp2.getResolutionMethod(),MassAnalyzer::FWHM);
	TEST_EQUAL(tmp2.getResolutionType(),MassAnalyzer::CONSTANT);
	TEST_EQUAL(tmp2.getScanDirection(),MassAnalyzer::UP);
	TEST_EQUAL(tmp2.getScanLaw(),MassAnalyzer::LINEAR);
	TEST_REAL_SIMILAR(tmp2.getScanRate(),47.15);
	TEST_REAL_SIMILAR(tmp2.getScanTime(),47.16);
	TEST_REAL_SIMILAR(tmp2.getTOFTotalPathLength(),47.17);
	TEST_EQUAL(tmp2.getTandemScanMethod(),MassAnalyzer::PRODUCTIONSCAN);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)

	tmp2 = MassAnalyzer();
	TEST_EQUAL(tmp2.getType(),MassAnalyzer::ANALYZERNULL);
	TEST_REAL_SIMILAR(tmp2.getAccuracy(),0.0);
	TEST_EQUAL(tmp2.getFinalMSExponent(),0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWidth(),0.0);
	TEST_REAL_SIMILAR(tmp2.getMagneticFieldStrength(),0.0);
	TEST_EQUAL(tmp2.getReflectronState(),MassAnalyzer::REFLSTATENULL);
	TEST_REAL_SIMILAR(tmp2.getResolution(),0.0);
	TEST_EQUAL(tmp2.getResolutionMethod(),MassAnalyzer::RESMETHNULL);
	TEST_EQUAL(tmp2.getResolutionType(),MassAnalyzer::RESTYPENULL);
	TEST_EQUAL(tmp2.getScanDirection(),MassAnalyzer::SCANDIRNULL);
	TEST_EQUAL(tmp2.getScanLaw(),MassAnalyzer::SCANLAWNULL);
	TEST_REAL_SIMILAR(tmp2.getScanRate(),0.0);
	TEST_REAL_SIMILAR(tmp2.getScanTime(),0.0);
	TEST_REAL_SIMILAR(tmp2.getTOFTotalPathLength(),0.0);
	TEST_EQUAL(tmp2.getTandemScanMethod(),MassAnalyzer::TANDEMNULL);
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getOrder(),0)
END_SECTION

START_SECTION(bool operator== (const MassAnalyzer& rhs) const)
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
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION(bool operator!= (const MassAnalyzer& rhs) const)
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
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit!=empty,true);
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



