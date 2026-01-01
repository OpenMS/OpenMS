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

#include <unordered_set>
#include <unordered_map>

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
  TEST_EQUAL(tmp.getType(),MassAnalyzer::ANALYZERNULL);
END_SECTION

START_SECTION((ReflectronState getReflectronState() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::REFLSTATENULL);
END_SECTION

START_SECTION((ResolutionMethod getResolutionMethod() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::RESMETHNULL);
END_SECTION

START_SECTION((ResolutionType getResolutionType() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::RESTYPENULL);
END_SECTION

START_SECTION((ScanDirection getScanDirection() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::SCANDIRNULL);
END_SECTION

START_SECTION((ScanLaw getScanLaw() const))
  MassAnalyzer tmp;
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::SCANLAWNULL);
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
  tmp.setType(MassAnalyzer::QUADRUPOLE);
  TEST_EQUAL(tmp.getType(),MassAnalyzer::QUADRUPOLE);
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
  tmp.setReflectronState(MassAnalyzer::ON);
  TEST_EQUAL(tmp.getReflectronState(),MassAnalyzer::ON);
END_SECTION

START_SECTION((void setResolution(double resolution)))
  MassAnalyzer tmp;
  tmp.setResolution(47.11);
  TEST_REAL_SIMILAR(tmp.getResolution(),47.11);
END_SECTION

START_SECTION((void setResolutionMethod(ResolutionMethod resolution_method)))
  MassAnalyzer tmp;
  tmp.setResolutionMethod(MassAnalyzer::FWHM);
  TEST_EQUAL(tmp.getResolutionMethod(),MassAnalyzer::FWHM);
END_SECTION

START_SECTION((void setResolutionType(ResolutionType resolution_type)))
  MassAnalyzer tmp;
  tmp.setResolutionType(MassAnalyzer::CONSTANT);
  TEST_EQUAL(tmp.getResolutionType(),MassAnalyzer::CONSTANT);
END_SECTION

START_SECTION((void setScanDirection(ScanDirection scan_direction)))
  MassAnalyzer tmp;
  tmp.setScanDirection(MassAnalyzer::UP);
  TEST_EQUAL(tmp.getScanDirection(),MassAnalyzer::UP);
END_SECTION

START_SECTION((void setScanLaw(ScanLaw scan_law)))
  MassAnalyzer tmp;
  tmp.setScanLaw(MassAnalyzer::LINEAR);
  TEST_EQUAL(tmp.getScanLaw(),MassAnalyzer::LINEAR);
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
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_EQUAL(tmp2.getOrder(),45)
END_SECTION

START_SECTION((MassAnalyzer& operator= (const MassAnalyzer& source)))
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
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_EQUAL(tmp2.getOrder(),0)
END_SECTION

START_SECTION((bool operator== (const MassAnalyzer& rhs) const))
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
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit==empty,false);
END_SECTION

START_SECTION((bool operator!= (const MassAnalyzer& rhs) const))
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
	
  edit = empty;
  edit.setOrder(45);
	TEST_EQUAL(edit!=empty,true);
END_SECTION

START_SECTION((static StringList getAllNamesOfAnalyzerType()))
  StringList names = MassAnalyzer::getAllNamesOfAnalyzerType();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_ANALYZERTYPE);
  TEST_EQUAL(names[MassAnalyzer::QUADRUPOLE], "Quadrupole");
  TEST_EQUAL(names[MassAnalyzer::ORBITRAP], "Orbitrap");
END_SECTION

START_SECTION((static StringList getAllNamesOfResolutionMethod()))
  StringList names = MassAnalyzer::getAllNamesOfResolutionMethod();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_RESOLUTIONMETHOD);
  TEST_EQUAL(names[MassAnalyzer::FWHM], "Full width at half max");
END_SECTION

START_SECTION((static StringList getAllNamesOfResolutionType()))
  StringList names = MassAnalyzer::getAllNamesOfResolutionType();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_RESOLUTIONTYPE);
  TEST_EQUAL(names[MassAnalyzer::CONSTANT], "Constant");
END_SECTION

START_SECTION((static StringList getAllNamesOfScanDirection()))
  StringList names = MassAnalyzer::getAllNamesOfScanDirection();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_SCANDIRECTION);
  TEST_EQUAL(names[MassAnalyzer::UP], "Up");
END_SECTION

START_SECTION((static StringList getAllNamesOfScanLaw()))
  StringList names = MassAnalyzer::getAllNamesOfScanLaw();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_SCANLAW);
  TEST_EQUAL(names[MassAnalyzer::LINEAR], "Linar");
END_SECTION

START_SECTION((static StringList getAllNamesOfReflectronState()))
  StringList names = MassAnalyzer::getAllNamesOfReflectronState();
  TEST_EQUAL(names.size(), MassAnalyzer::SIZE_OF_REFLECTRONSTATE);
  TEST_EQUAL(names[MassAnalyzer::ON], "On");
END_SECTION

START_SECTION(([EXTRA] std::hash<MassAnalyzer>))
{
  // Test that equal objects have equal hashes
  MassAnalyzer ma1, ma2;
  std::hash<MassAnalyzer> hasher;

  // Default constructed objects should have equal hashes
  TEST_EQUAL(hasher(ma1), hasher(ma2))

  // Set all fields to the same values
  ma1.setType(MassAnalyzer::QUADRUPOLE);
  ma1.setResolutionMethod(MassAnalyzer::FWHM);
  ma1.setResolutionType(MassAnalyzer::CONSTANT);
  ma1.setScanDirection(MassAnalyzer::UP);
  ma1.setScanLaw(MassAnalyzer::LINEAR);
  ma1.setReflectronState(MassAnalyzer::ON);
  ma1.setResolution(47.14);
  ma1.setAccuracy(47.11);
  ma1.setScanRate(47.15);
  ma1.setScanTime(47.16);
  ma1.setTOFTotalPathLength(47.17);
  ma1.setIsolationWidth(47.12);
  ma1.setMagneticFieldStrength(47.13);
  ma1.setFinalMSExponent(47);
  ma1.setOrder(45);

  ma2.setType(MassAnalyzer::QUADRUPOLE);
  ma2.setResolutionMethod(MassAnalyzer::FWHM);
  ma2.setResolutionType(MassAnalyzer::CONSTANT);
  ma2.setScanDirection(MassAnalyzer::UP);
  ma2.setScanLaw(MassAnalyzer::LINEAR);
  ma2.setReflectronState(MassAnalyzer::ON);
  ma2.setResolution(47.14);
  ma2.setAccuracy(47.11);
  ma2.setScanRate(47.15);
  ma2.setScanTime(47.16);
  ma2.setTOFTotalPathLength(47.17);
  ma2.setIsolationWidth(47.12);
  ma2.setMagneticFieldStrength(47.13);
  ma2.setFinalMSExponent(47);
  ma2.setOrder(45);

  // Equal objects should have equal hashes
  TEST_EQUAL(ma1 == ma2, true)
  TEST_EQUAL(hasher(ma1), hasher(ma2))

  // Test that different objects have different hashes (not guaranteed but highly likely)
  MassAnalyzer ma3;
  ma3.setType(MassAnalyzer::TOF);
  TEST_NOT_EQUAL(hasher(ma1), hasher(ma3))

  // Test use in unordered_set
  std::unordered_set<MassAnalyzer> ma_set;
  ma_set.insert(ma1);
  ma_set.insert(ma2); // Should not add a duplicate
  ma_set.insert(ma3);
  TEST_EQUAL(ma_set.size(), 2) // ma1 and ma2 are equal, so only 2 elements

  // Test use in unordered_map
  std::unordered_map<MassAnalyzer, int> ma_map;
  ma_map[ma1] = 1;
  ma_map[ma2] = 2; // Should overwrite ma1's value since they're equal
  ma_map[ma3] = 3;
  TEST_EQUAL(ma_map.size(), 2)
  TEST_EQUAL(ma_map[ma1], 2) // Value should be 2 since ma2 overwrote ma1

  // Test hash consistency (same object should produce same hash)
  std::size_t hash1 = hasher(ma1);
  std::size_t hash2 = hasher(ma1);
  TEST_EQUAL(hash1, hash2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



