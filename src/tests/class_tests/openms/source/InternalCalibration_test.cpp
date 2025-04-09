// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

START_TEST(InternalCalibration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

InternalCalibration* ptr = nullptr;
InternalCalibration* nullPointer = nullptr;
START_SECTION(InternalCalibration())
{
  ptr = new InternalCalibration();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~InternalCalibration())
{
  delete ptr;
}
END_SECTION


START_SECTION(Size fillCalibrants(const PeakMap exp, const std::vector<InternalCalibration::LockMass>& ref_masses, double tol_ppm, bool lock_require_mono, bool lock_require_iso, CalibrationData& failed_lock_masses, bool verbose = true))
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("InternalCalibration_2_lockmass.mzML.gz"), exp);
  std::vector<InternalCalibration::LockMass> ref_masses;
  
  ref_masses.push_back(InternalCalibration::LockMass(327.25353, 1, 1));
  ref_masses.push_back(InternalCalibration::LockMass(362.29065, 1, 1));
  ref_masses.push_back(InternalCalibration::LockMass(680.48022, 1, 1));

  InternalCalibration ic;
  CalibrationData failed_locks;
  Size cal_count = ic.fillCalibrants(exp, ref_masses, 25.0, true, false, failed_locks, true); // no 'require_iso', since the example data has really high C13 mass error (up to 7ppm to +1 iso)

  TEST_EQUAL(cal_count, 21 * 3); // 21 MS1 scans, 3 calibrants each

END_SECTION

std::vector<PeptideIdentification> peps;
std::vector<ProteinIdentification> prots;
IdXMLFile().load(File::find("./examples/BSA/BSA1_OMSSA.idXML"), prots, peps);

START_SECTION(Size fillCalibrants(const FeatureMap& fm, double tol_ppm))
  FeatureMap fm;
  fm.setUnassignedPeptideIdentifications(peps);

  InternalCalibration ic;
  Size cal_count = ic.fillCalibrants(fm, 100.0);
  TEST_EQUAL(cal_count, 44); // all pep IDs

  cal_count = ic.fillCalibrants(fm, 10.0);
  TEST_EQUAL(cal_count, 37);  // a few outliers IDs removed

END_SECTION

START_SECTION(Size fillCalibrants(const std::vector<PeptideIdentification>& pep_ids, double tol_ppm))
  InternalCalibration ic;
  Size cal_count = ic.fillCalibrants(peps, 100.0);
  TEST_EQUAL(cal_count, 44);

  cal_count = ic.fillCalibrants(peps, 10.0);
  TEST_EQUAL(cal_count, 37);

  TEST_EQUAL(ic.getCalibrationPoints().size(), cal_count)

END_SECTION

START_SECTION(const CalibrationData& getCalibrationPoints() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(bool calibrate(PeakMap& exp, const IntList& target_mslvl, MZTrafoModel::MODELTYPE model_type, double rt_chunk, bool use_RANSAC, double post_ppm_median, double post_ppm_MAD, const String& file_models, const String& file_residuals))
  InternalCalibration ic;
  ic.fillCalibrants(peps, 3.0);
  PeakMap exp;
  MzMLFile().load(File::find("./examples/BSA/BSA1.mzML"), exp);
  MZTrafoModel::setRANSACSeed(0);
  MZTrafoModel::setRANSACParams(Math::RANSACParam(2, 1000, 1.0, 30, true));
  bool success = ic.calibrate(exp, std::vector<Int>(1, 1), MZTrafoModel::LINEAR, -1, true, 1.0, 1.0);
  TEST_EQUAL(success, true)
END_SECTION

PeakMap::SpectrumType spec;
spec.push_back(Peak1D(250.0, 1000.0));
spec.push_back(Peak1D(500.0, 1000.0));
spec.push_back(Peak1D(750.0, 1000.0));
spec.push_back(Peak1D(1000.0, 1000.0));
std::vector<Precursor> pcs;
Precursor pc;
pc.setMZ(123.0);
pcs.push_back(pc);
pc.setMZ(456.0);
pcs.push_back(pc);
spec.setPrecursors(pcs);

START_SECTION(static void applyTransformation(std::vector<Precursor>& pcs, const MZTrafoModel& trafo))
  MZTrafoModel trafo;
  trafo.setCoefficients(-100.0, 0.0, 0.0);
  std::vector<Precursor> pcs2 = pcs;
  InternalCalibration::applyTransformation(pcs2, trafo);
  TEST_REAL_SIMILAR(pcs2[0].getMZ(), pcs[0].getMZ() - Math::ppmToMass(-100.0, 123.0));
  TEST_REAL_SIMILAR(pcs2[1].getMZ(), pcs[1].getMZ() - Math::ppmToMass(-100.0, 456.0));
  //test for meta value "mz_raw"
  ABORT_IF(!pcs2[0].metaValueExists("mz_raw"));
  TEST_REAL_SIMILAR(pcs2[0].getMetaValue("mz_raw"), 123);
  ABORT_IF(!pcs2[1].metaValueExists("mz_raw"));
  TEST_REAL_SIMILAR(pcs2[1].getMetaValue("mz_raw"), 456);

END_SECTION

START_SECTION(static void applyTransformation(PeakMap::SpectrumType& spec, const IntList& target_mslvl, const MZTrafoModel& trafo))
  MZTrafoModel trafo;
  trafo.setCoefficients(-100.0, 0.0, 0.0);
  PeakMap::SpectrumType spec2 = spec;
  TEST_EQUAL(spec, spec2);
  InternalCalibration::applyTransformation(spec2, std::vector<Int>(1, 1), trafo);
  TEST_NOT_EQUAL(spec, spec2);
  TEST_REAL_SIMILAR(spec2[0].getMZ(), spec[0].getMZ() - Math::ppmToMass(-100.0, 250.0));
  TEST_REAL_SIMILAR(spec2[1].getMZ(), spec[1].getMZ() - Math::ppmToMass(-100.0, 500.0));
  TEST_EQUAL(spec2.getPrecursors()[0], pcs[0]); // unchanged, since PCs belong to MS-level 0
  TEST_EQUAL(spec2.getPrecursors()[1], pcs[1]); // unchanged, since PCs belong to MS-level 0

  spec2 = spec;
  spec2.setMSLevel(2);
  PeakMap::SpectrumType spec2_noPC = spec2;
  spec2_noPC.getPrecursors().resize(0); // remove PC's
  InternalCalibration::applyTransformation(spec2, std::vector<Int>(1, 1), trafo);
  TEST_REAL_SIMILAR(spec2.getPrecursors()[0].getMZ(), pcs[0].getMZ() - Math::ppmToMass(-100.0, 123.0));
  TEST_REAL_SIMILAR(spec2.getPrecursors()[1].getMZ(), pcs[1].getMZ() - Math::ppmToMass(-100.0, 456.0));
  spec2.getPrecursors().resize(0); // remove PC's
  TEST_EQUAL(spec2_noPC, spec2); // everything else should be unchanged

END_SECTION

START_SECTION(static void applyTransformation(PeakMap& exp, const IntList& target_mslvl, const MZTrafoModel& trafo))
  MZTrafoModel trafo;
  trafo.setCoefficients(-100.0, 0.0, 0.0); // observed m/z are 100ppm lower than reference
  PeakMap::SpectrumType spec2 = spec;
  spec2.setMSLevel(2);      // will not be calibrated, except for its PC
  PeakMap exp;
  exp.addSpectrum(spec);
  exp.addSpectrum(spec2);
  exp.addSpectrum(spec);
  
  InternalCalibration::applyTransformation(exp, std::vector<Int>(1, 1), trafo);
  TEST_NOT_EQUAL(exp[0], spec);
  TEST_REAL_SIMILAR(exp[0][0].getMZ(), spec[0].getMZ() + Math::ppmToMass(-1 * -100.0, 250.0));
  TEST_REAL_SIMILAR(exp[0][1].getMZ(), spec[1].getMZ() + Math::ppmToMass(-1 *-100.0, 500.0));
  TEST_REAL_SIMILAR(spec.getPrecursors()[0].getMZ(), exp[0].getPrecursors()[0].getMZ());
  TEST_REAL_SIMILAR(spec.getPrecursors()[1].getMZ(), exp[0].getPrecursors()[1].getMZ());

  TEST_NOT_EQUAL(exp[1], spec2);
  TEST_REAL_SIMILAR(exp[1][0].getMZ(), spec2[0].getMZ());
  TEST_REAL_SIMILAR(exp[1][1].getMZ(), spec2[1].getMZ());
  TEST_REAL_SIMILAR(spec2.getPrecursors()[0].getMZ(), exp[1].getPrecursors()[0].getMZ() + Math::ppmToMass(-100.0, 123.0));
  TEST_REAL_SIMILAR(spec2.getPrecursors()[1].getMZ(), exp[1].getPrecursors()[1].getMZ() + Math::ppmToMass(-100.0, 456.0));
  
  TEST_NOT_EQUAL(exp[2], spec);
  TEST_REAL_SIMILAR(exp[2][0].getMZ(), spec[0].getMZ() + Math::ppmToMass(-1 *-100.0, 250.0));
  TEST_REAL_SIMILAR(exp[2][1].getMZ(), spec[1].getMZ() + Math::ppmToMass(-1 *-100.0, 500.0));
  TEST_REAL_SIMILAR(spec.getPrecursors()[0].getMZ(), exp[2].getPrecursors()[0].getMZ());
  TEST_REAL_SIMILAR(spec.getPrecursors()[1].getMZ(), exp[2].getPrecursors()[1].getMZ());

END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


