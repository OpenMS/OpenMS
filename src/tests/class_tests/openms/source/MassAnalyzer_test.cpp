// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/MassAnalyzer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassAnalyzer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassAnalyzer* ptr = nullptr;
MassAnalyzer* nullPointer = nullptr;
START_SECTION((MassAnalyzer()))
	ptr = new MassAnalyzer();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MassAnalyzer()))
	delete ptr;
END_SECTION

START_SECTION((Int getOrder() const))
	MassAnalyzer tmp;
	TEST_EQUAL(tmp.getOrder(),0)
END_SECTION

START_SECTION((void setOrder(Int order)))
	MassAnalyzer tmp;
	tmp.setOrder(4711);
	TEST_EQUAL(tmp.getOrder(),4711)
END_SECTION

START_SECTION((AnalyzerType getType() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getType(),MassAnalyzer::AnalyzerType::ANALYZERNULL);
END_SECTION

START_SECTION((ReflectronState getReflectronState() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ReflectronState::REFLSTATENULL);
END_SECTION

START_SECTION((ResolutionMethod getResolutionMethod() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::ResolutionMethod::RESMETHNULL);
END_SECTION

START_SECTION((ResolutionType getResolutionType() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::ResolutionType::RESTYPENULL);
END_SECTION

START_SECTION((ScanDirection getScanDirection() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::ScanDirection::SCANDIRNULL);
END_SECTION

START_SECTION((ScanLaw getScanLaw() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::ScanLaw::SCANLAWNULL);
END_SECTION

START_SECTION((Int getFinalMSExponent() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getFinalMSExponent(),0);
END_SECTION

START_SECTION((double getAccuracy() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getAccuracy(),0.0);
END_SECTION

START_SECTION((double getIsolationWidth() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWidth(),0.0);
END_SECTION

START_SECTION((double getMagneticFieldStrength() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getMagneticFieldStrength(),0.0);
END_SECTION

START_SECTION((double getResolution() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getResolution(),0.0);
END_SECTION

START_SECTION((double getScanRate() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getScanRate(),0.0);
END_SECTION

START_SECTION((double getScanTime() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getScanTime(),0.0);
END_SECTION

START_SECTION((double getTOFTotalPathLength() const ))
  MassAnalyzer tmp;
  TEST_REAL_SIMILAR(tmp.getTOFTotalPathLength(),0.0);
END_SECTION

START_SECTION((void setType(AnalyzerType type)))
  MassAnalyzer tmp;
  tmp.setType(MassAnalyzer::AnalyzerType::QUADRUPOLE);
  TEST_EQUAL(tmp.getType(),MassAnalyzer::AnalyzerType::QUADRUPOLE);
END_SECTION

START_SECTION((void setAccuracy(double accuracy)))
  MassAnalyzer tmp;
  tmp.setAccuracy(47.11);
  TEST_REAL_SIMILAR(tmp.getAccuracy(),47.11);
END_SECTION

START_SECTION((void setFinalMSExponent(Int final_MS_exponent)))
  MassAnalyzer tmp;
  tmp.setFinalMSExponent(47);
  TEST_EQUAL(tmp.getFinalMSExponent(),47);
END_SECTION

START_SECTION((void setIsolationWidth(double isolation_width)))
  MassAnalyzer tmp;
  tmp.setIsolationWidth(47.11);
  TEST_REAL_SIMILAR(tmp.getIsolationWidth(),47.11);
END_SECTION

START_SECTION((void setMagneticFieldStrength(double magnetic_field_strength)))
  MassAnalyzer tmp;
  tmp.setMagneticFieldStrength(47.11);
  TEST_REAL_SIMILAR(tmp.getMagneticFieldStrength(),47.11);
END_SECTION

START_SECTION((void setReflectronState(ReflectronState reflecton_state)))
  MassAnalyzer tmp;
  tmp.setReflectronState(MassAnalyzer::ReflectronState::ON);
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ReflectronState::ON);
END_SECTION

START_SECTION((void setResolution(double resolution)))
  MassAnalyzer tmp;
  tmp.setResolution(47.11);
  TEST_REAL_SIMILAR(tmp.getResolution(),47.11);
END_SECTION

START_SECTION((void setResolutionMethod(ResolutionMethod resolution_method)))
  MassAnalyzer tmp;
  tmp.setResolutionMethod(MassAnalyzer::ResolutionMethod::FWHM);
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::ResolutionMethod::FWHM);
END_SECTION

START_SECTION((void setResolutionType(ResolutionType resolution_type)))
  MassAnalyzer tmp;
  tmp.setResolutionType(MassAnalyzer::ResolutionType::CONSTANT);
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::ResolutionType::CONSTANT);
END_SECTION

START_SECTION((void setScanDirection(ScanDirection scan_direction)))
  MassAnalyzer tmp;
  tmp.setScanDirection(MassAnalyzer::ScanDirection::UP);
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::ScanDirection::UP);
END_SECTION

START_SECTION((void setScanLaw(ScanLaw scan_law)))
  MassAnalyzer tmp;
  tmp.setScanLaw(MassAnalyzer::ScanLaw::LINEAR);
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::ScanLaw::LINEAR);
END_SECTION

START_SECTION((void setScanRate(double scan_rate)))
  MassAnalyzer tmp;
  tmp.setScanRate(47.11);
  TEST_REAL_SIMILAR(tmp.getScanRate(),47.11);
END_SECTION

START_SECTION((void setScanTime(double scan_time)))
  MassAnalyzer tmp;
  tmp.setScanTime(47.11);
  TEST_REAL_SIMILAR(tmp.getScanTime(),47.11);
END_SECTION

START_SECTION((void setTOFTotalPathLength(double TOF_total_path_length)))
  MassAnalyzer tmp;
  tmp.setTOFTotalPathLength(47.11);
  TEST_REAL_SIMILAR(tmp.getTOFTotalPathLength(),47.11);
END_SECTION

START_SECTION((MassAnalyzer(const MassAnalyzer& source)))
  MassAnalyzer tmp;
  tmp.setType(MassAnalyzer::AnalyzerType::QUADRUPOLE);
  tmp.setAccuracy(47.11);
  tmp.setFinalMSExponent(47);
  tmp.setIsolationWidth(47.12);
  tmp.setMagneticFieldStrength(47.13);
  tmp.setReflectronState(MassAnalyzer::ReflectronState::ON);
  tmp.setResolution(47.14);
  tmp.setResolutionMethod(MassAnalyzer::ResolutionMethod::FWHM);
  tmp.setResolutionType(MassAnalyzer::ResolutionType::CONSTANT);
  tmp.setScanDirection(MassAnalyzer::ScanDirection::UP);
  tmp.setScanLaw(MassAnalyzer::ScanLaw::LINEAR);
  tmp.setScanRate(47.15);
  tmp.setScanTime(47.16);
  tmp.setTOFTotalPathLength(47.17);
	tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
	
	MassAnalyzer tmp2(tmp);
	TEST_EQUAL(tmp.getType(),MassAnalyzer::AnalyzerType::QUADRUPOLE);
	TEST_REAL_SIMILAR(tmp.getAccuracy(),47.11);
	TEST_EQUAL(tmp.getFinalMSExponent(),47);
	TEST_REAL_SIMILAR(tmp.getIsolationWidth(),47.12);
	TEST_REAL_SIMILAR(tmp.getMagneticFieldStrength(),47.13);
	TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ReflectronState::ON);
	TEST_REAL_SIMILAR(tmp.getResolution(),47.14);
	TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::ResolutionMethod::FWHM);
	TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::ResolutionType::CONSTANT);
	TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::ScanDirection::UP);
	TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::ScanLaw::LINEAR);
	TEST_REAL_SIMILAR(tmp.getScanRate(),47.15);
	TEST_REAL_SIMILAR(tmp.getScanTime(),47.16);
	TEST_REAL_SIMILAR(tmp.getTOFTotalPathLength(),47.17);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION((MassAnalyzer& operator= (const MassAnalyzer& source)))
  MassAnalyzer tmp;
  
  tmp.setType(MassAnalyzer::AnalyzerType::QUADRUPOLE);
  tmp.setAccuracy(47.11);
  tmp.setFinalMSExponent(47);
  tmp.setIsolationWidth(47.12);
  tmp.setMagneticFieldStrength(47.13);
  tmp.setReflectronState(MassAnalyzer::ReflectronState::ON);
  tmp.setResolution(47.14);
  tmp.setResolutionMethod(MassAnalyzer::ResolutionMethod::FWHM);
  tmp.setResolutionType(MassAnalyzer::ResolutionType::CONSTANT);
  tmp.setScanDirection(MassAnalyzer::ScanDirection::UP);
  tmp.setScanLaw(MassAnalyzer::ScanLaw::LINEAR);
  tmp.setScanRate(47.15);
  tmp.setScanTime(47.16);
  tmp.setTOFTotalPathLength(47.17);
	tmp.setMetaValue("label",String("label"));
  tmp.setOrder(45);
	
	MassAnalyzer tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getType(),MassAnalyzer::AnalyzerType::QUADRUPOLE);
	TEST_REAL_SIMILAR(tmp2.getAccuracy(),47.11);
	TEST_EQUAL(tmp2.getFinalMSExponent(),47);
	TEST_REAL_SIMILAR(tmp2.getIsolationWidth(),47.12);
	TEST_REAL_SIMILAR(tmp2.getMagneticFieldStrength(),47.13);
	TEST_EQUAL(tmp2.getReflectronState(),MassAnalyzer::ReflectronState::ON);
	TEST_REAL_SIMILAR(tmp2.getResolution(),47.14);
	TEST_EQUAL(tmp2.getResolutionMethod(),MassAnalyzer::ResolutionMethod::FWHM);
	TEST_EQUAL(tmp2.getResolutionType(),MassAnalyzer::ResolutionType::CONSTANT);
	TEST_EQUAL(tmp2.getScanDirection(),MassAnalyzer::ScanDirection::UP);
	TEST_EQUAL(tmp2.getScanLaw(),MassAnalyzer::ScanLaw::LINEAR);
	TEST_REAL_SIMILAR(tmp2.getScanRate(),47.15);
	TEST_REAL_SIMILAR(tmp2.getScanTime(),47.16);
	TEST_REAL_SIMILAR(tmp2.getTOFTotalPathLength(),47.17);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)

	tmp2 = MassAnalyzer();
	TEST_EQUAL(tmp2.getType(),MassAnalyzer::AnalyzerType::ANALYZERNULL);
	TEST_REAL_SIMILAR(tmp2.getAccuracy(),0.0);
	TEST_EQUAL(tmp2.getFinalMSExponent(),0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWidth(),0.0);
	TEST_REAL_SIMILAR(tmp2.getMagneticFieldStrength(),0.0);
	TEST_EQUAL(tmp2.getReflectronState(),MassAnalyzer::ReflectronState::REFLSTATENULL);
	TEST_REAL_SIMILAR(tmp2.getResolution(),0.0);
	TEST_EQUAL(tmp2.getResolutionMethod(),MassAnalyzer::ResolutionMethod::RESMETHNULL);
	TEST_EQUAL(tmp2.getResolutionType(),MassAnalyzer::ResolutionType::RESTYPENULL);
	TEST_EQUAL(tmp2.getScanDirection(),MassAnalyzer::ScanDirection::SCANDIRNULL);
	TEST_EQUAL(tmp2.getScanLaw(),MassAnalyzer::ScanLaw::SCANLAWNULL);
	TEST_REAL_SIMILAR(tmp2.getScanRate(),0.0);
	TEST_REAL_SIMILAR(tmp2.getScanTime(),0.0);
	TEST_REAL_SIMILAR(tmp2.getTOFTotalPathLength(),0.0);
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getOrder(),0)
END_SECTION

START_SECTION((bool operator== (const MassAnalyzer& rhs) const))
  MassAnalyzer edit, empty;
  
	TEST_EQUAL(edit==empty,true);

	edit=empty;
	edit.setType(MassAnalyzer::AnalyzerType::QUADRUPOLE);
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
	edit.setReflectronState(MassAnalyzer::ReflectronState::ON);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setResolution(47.14);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setResolutionMethod(MassAnalyzer::ResolutionMethod::FWHM);
	TEST_EQUAL(edit==empty,false);

	edit=empty;
	edit.setResolutionType(MassAnalyzer::ResolutionType::CONSTANT);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanDirection(MassAnalyzer::ScanDirection::UP);
	TEST_EQUAL(edit==empty,false);
	
	edit=empty;
	edit.setScanLaw(MassAnalyzer::ScanLaw::LINEAR);
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
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!= (const MassAnalyzer& rhs) const))
  MassAnalyzer edit, empty;
  
	TEST_EQUAL(edit!=empty,false);

	edit=empty;
	edit.setType(MassAnalyzer::AnalyzerType::QUADRUPOLE);
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
	edit.setReflectronState(MassAnalyzer::ReflectronState::ON);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setResolution(47.14);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setResolutionMethod(MassAnalyzer::ResolutionMethod::FWHM);
	TEST_EQUAL(edit!=empty,true);

	edit=empty;
	edit.setResolutionType(MassAnalyzer::ResolutionType::CONSTANT);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanDirection(MassAnalyzer::ScanDirection::UP);
	TEST_EQUAL(edit!=empty,true);
	
	edit=empty;
	edit.setScanLaw(MassAnalyzer::ScanLaw::LINEAR);
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
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit!=empty,true);
END_SECTION

START_SECTION((static StringList getAllNamesOfAnalyzerType()))
  StringList names = MassAnalyzer::getAllNamesOfAnalyzerType();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_ANALYZERTYPE);
  TEST_EQUAL(names[MassAnalyzer::AnalyzerType::QUADRUPOLE], "Quadrupole");
  TEST_EQUAL(names[MassAnalyzer::AnalyzerType::ORBITRAP], "Orbitrap");
END_SECTION

START_SECTION((static StringList getAllNamesOfResolutionMethod()))
  StringList names = MassAnalyzer::getAllNamesOfResolutionMethod();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
  TEST_EQUAL(names[MassAnalyzer::ResolutionMethod::FWHM], "Full width at half max");
END_SECTION

START_SECTION((static StringList getAllNamesOfResolutionType()))
  StringList names = MassAnalyzer::getAllNamesOfResolutionType();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_RESOLUTIONTYPE);
  TEST_EQUAL(names[MassAnalyzer::ResolutionType::CONSTANT], "Constant");
END_SECTION

START_SECTION((static StringList getAllNamesOfScanDirection()))
  StringList names = MassAnalyzer::getAllNamesOfScanDirection();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_SCANDIRECTION);
  TEST_EQUAL(names[MassAnalyzer::ScanDirection::UP], "Up");
END_SECTION

START_SECTION((static StringList getAllNamesOfScanLaw()))
  StringList names = MassAnalyzer::getAllNamesOfScanLaw();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_SCANLAW);
  TEST_EQUAL(names[MassAnalyzer::ScanLaw::LINEAR], "Linar");
END_SECTION

START_SECTION((static StringList getAllNamesOfReflectronState()))
  StringList names = MassAnalyzer::getAllNamesOfReflectronState();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_REFLECTRONSTATE);
  TEST_EQUAL(names[MassAnalyzer::ReflectronState::ON], "On");
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



