// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QFile>

using namespace OpenMS;
using namespace std;

///////////////////////////


void cmpDataIntensity(const MSExperiment& exp1, const MSExperiment& exp2, double abs_tol = 1e-5, double rel_tol = 1+1e-5)
{
  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  TOLERANCE_ABSOLUTE(abs_tol)
  TOLERANCE_RELATIVE(rel_tol)
  for (Size i = 0; i < exp1.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp1.getSpectra()[i].size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp1.getSpectra()[i].size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      auto a = exp1.getSpectra()[i][k].getIntensity();
      auto b = exp2.getSpectra()[i][k].getIntensity();
      // avoid the console clutter if nothing interesing happens
      if (!TEST::isRealSimilar(a, b)) TEST_REAL_SIMILAR(a, b)
    }
  }

  for (Size i = 0; i < exp1.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp1.getChromatograms()[i].size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp1.getChromatograms()[i].size(); k++)
    {
      auto a = exp1.getChromatograms()[i][k].getIntensity();
      auto b = exp2.getChromatograms()[i][k].getIntensity();
      // avoid the console clutter if nothing interesing happens
      if (!TEST::isRealSimilar(a, b)) TEST_REAL_SIMILAR(a, b)
    }
  }
}

void cmpDataMZ(const MSExperiment& exp1, const MSExperiment& exp2, double abs_tol = 1e-5, double rel_tol = 1+1e-5)
{
  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  TOLERANCE_ABSOLUTE(abs_tol)
  TOLERANCE_RELATIVE(rel_tol)
  for (Size i = 0; i < exp1.getNrSpectra(); i++)
  {
    TEST_EQUAL(exp1.getSpectra()[i].size(), exp2.getSpectra()[i].size())
    for (Size k = 0; k < exp1.getSpectra()[i].size(); k++)
    {
      // slof is no good for values smaller than 5
      // if (exp.getSpectrum(i)[k].getIntensity() < 1.0) {continue;} 
      auto a = exp1.getSpectra()[i][k].getMZ();
      auto b = exp2.getSpectra()[i][k].getMZ();
      // avoid the console clutter if nothing interesing happens
      if (!TEST::isRealSimilar(a, b)) TEST_REAL_SIMILAR(a, b)
    }
  }
}

void cmpDataRT(const MSExperiment& exp1, const MSExperiment& exp2, double abs_tol = 1e-5, double rel_tol = 1+1e-5)
{
  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  TOLERANCE_ABSOLUTE(abs_tol)
  TOLERANCE_RELATIVE(rel_tol)
  for (Size i = 0; i < exp1.getNrChromatograms(); i++)
  {
    TEST_EQUAL(exp1.getChromatograms()[i].size() == exp2.getChromatograms()[i].size(), true)
    for (Size k = 0; k < exp1.getChromatograms()[i].size(); k++)
    {
      auto a = exp1.getChromatograms()[i][k].getRT();
      auto b = exp2.getChromatograms()[i][k].getRT();
      // avoid the console clutter if nothing interesing happens
      if (!TEST::isRealSimilar(a, b)) TEST_REAL_SIMILAR(a, b)
    }
  }
}


START_TEST(SqMassFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SqMassFile* ptr = nullptr;
SqMassFile* nullPointer = nullptr;
START_SECTION((SqMassFile()))
  ptr = new SqMassFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~SqMassFile()))
  delete ptr;
END_SECTION

TOLERANCE_RELATIVE(1.0005)

START_SECTION(void load(const String& filename, MapType& map))
{
  MSExperiment exp;
  SqMassFile().load(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), exp);

  MSExperiment exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLSqliteHandler_1.mzML"),exp2);

  TEST_EQUAL(exp.getNrSpectra(), exp2.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp2.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp2.getSpectra()[0], false) // no exact duplicate

  cout.precision(17);

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.

  cmpDataIntensity(exp, exp2, 1e-4, 1.001);
  cmpDataMZ(exp, exp2, 1e-5, 1.000001);  // less than 1ppm error for m/z 
  cmpDataRT(exp, exp2, 0.05, 1.000001);  // max 0.05 seconds error in RT

  // mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
  TEST_EQUAL(exp.getSqlRunID(), exp2.getSqlRunID())

}
END_SECTION

// reset error tolerances to default values
TOLERANCE_ABSOLUTE(1e-5)
TOLERANCE_RELATIVE(1+1e-5)

START_SECTION(void store(const String& filename, MapType& map))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLSqliteHandler_1.mzML"), exp_orig);

  SqMassFile::SqMassConfig config;
  config.use_lossy_numpress = false;
  config.linear_fp_mass_acc = -1;
  config.write_full_meta = false;

  SqMassFile file;
  file.setConfig(config);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename, exp_orig);

  MSExperiment exp;
  file.load(tmp_filename, exp);

  TEST_EQUAL(exp.getNrSpectra(), exp_orig.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp_orig.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp_orig.getSpectra()[0], false) // no exact duplicate

  MSExperiment exp2 = exp_orig;

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.
  TOLERANCE_ABSOLUTE(1e-4)
  TOLERANCE_RELATIVE(1.001) // 0.1 % error for intensity 

  // reset error tolerances to default values
  TOLERANCE_ABSOLUTE(1e-5)
  TOLERANCE_RELATIVE(1+1e-5)

  // since we specified no lossy compression, we expect high accuracy
  TOLERANCE_ABSOLUTE(1e-8)
  TOLERANCE_RELATIVE(1.00000001) 

  cmpDataIntensity(exp, exp2, 1e-8, 1.00000001);
  cmpDataMZ(exp, exp2, 1e-8, 1.00000001); 
  cmpDataRT(exp, exp2, 1e-8, 1.00000001); 

  // no 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, false)
}
END_SECTION

START_SECTION([EXTRA_LOSSY] void store(const String& filename, MapType& map))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLSqliteHandler_1.mzML"), exp_orig);

  SqMassFile::SqMassConfig config;
  config.use_lossy_numpress = true;
  config.linear_fp_mass_acc = 0.0001;
  config.write_full_meta = false;

  {
    Precursor p;
    std::set<Precursor::ActivationMethod> tmp = {Precursor::ActivationMethod::BIRD};
    p.setActivationMethods(tmp);
    p.setActivationEnergy(500);
    p.setCharge(4);
    p.setMZ(600);
    p.setIsolationWindowUpperOffset(7);
    p.setIsolationWindowLowerOffset(14);
    p.setDriftTime(0.5);
    p.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
    p.setMetaValue("peptide_sequence", "PEPTIDEK");

    std::vector<Precursor> prec;
    prec.push_back(p);
    exp_orig.getSpectrum(0).setPrecursors(prec);
    exp_orig.getChromatogram(0).setPrecursor(p);

    Product pr;
    // pr.setCharge(6);
    pr.setMZ(300);
    pr.setIsolationWindowUpperOffset(10);
    pr.setIsolationWindowLowerOffset(15);
    std::vector<Product> prod;
    prod.push_back(pr);
    exp_orig.getSpectrum(0).setProducts(prod);
    exp_orig.getChromatogram(0).setProduct(pr);
    TEST_REAL_SIMILAR(exp_orig.getSpectrum(0).getPrecursors()[0].getActivationEnergy(), 500.0);
  }

  SqMassFile file;
  file.setConfig(config);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename, exp_orig);

  MSExperiment exp;
  file.load(tmp_filename, exp);

  TEST_EQUAL(exp.getNrSpectra(), exp_orig.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp_orig.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp_orig.getSpectra()[0], false) // no exact duplicate

  {
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getRT(), exp_orig.getSpectrum(0).getRT())
    TEST_EQUAL(exp.getSpectrum(0).getNativeID(), exp_orig.getSpectrum(0).getNativeID())
    TEST_EQUAL(exp.getSpectrum(0).getMSLevel(), exp_orig.getSpectrum(0).getMSLevel())
    TEST_EQUAL(exp.getSpectrum(0).getInstrumentSettings().getPolarity(), exp_orig.getSpectrum(0).getInstrumentSettings().getPolarity())
    TEST_EQUAL(exp.getSpectrum(0).getProducts().size(), 1)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getProducts()[0].getMZ(), 300)
    // TEST_EQUAL(exp.getSpectrum(0).getProducts()[0].getCharge(), 6)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getProducts()[0].getIsolationWindowUpperOffset(), 10)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getProducts()[0].getIsolationWindowLowerOffset(), 15)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors().size(), 1)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getActivationEnergy(), 500)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors()[0].getCharge(), 4)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getDriftTime(), 0.5)
    // TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getDriftTimeUnit(), 0.5) // not implemented
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getMZ(), 600)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getIsolationWindowUpperOffset(), 7)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getIsolationWindowLowerOffset(), 14)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors()[0].metaValueExists("peptide_sequence"), true)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors()[0].getMetaValue("peptide_sequence"), "PEPTIDEK")
  }
  {
    TEST_EQUAL(exp.getChromatogram(0).getNativeID(), exp_orig.getChromatogram(0).getNativeID())
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getProduct().getMZ(), 300)
    // TEST_EQUAL(exp.getChromatogram(0).getProduct().getCharge(), 6)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getProduct().getIsolationWindowUpperOffset(), 10)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getProduct().getIsolationWindowLowerOffset(), 15)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getActivationEnergy(), 500)
    TEST_EQUAL(exp.getChromatogram(0).getPrecursor().getCharge(), 4)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getDriftTime(), 0.5)
    // TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getDriftTimeUnit(), 0.5) // not implemented
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getMZ(), 600)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getIsolationWindowUpperOffset(), 7)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getIsolationWindowLowerOffset(), 14)
    TEST_EQUAL(exp.getChromatogram(0).getPrecursor().metaValueExists("peptide_sequence"), true)
    TEST_EQUAL(exp.getChromatogram(0).getPrecursor().getMetaValue("peptide_sequence"), "PEPTIDEK")
  }

  MSExperiment exp2 = exp_orig;

  // should not give 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, false)

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.
  TOLERANCE_ABSOLUTE(1e-4)
  TOLERANCE_RELATIVE(1.001) // 0.1 % error for intensity 

  cmpDataIntensity(exp, exp2, 1e-4, 1.001);   // 0.1 % error for intensity 
  cmpDataMZ(exp, exp2, 1e-5, 1.000001);  // less than 1ppm error for m/z 
  cmpDataRT(exp, exp2, 0.05, 1.00000001);  // max 0.05 seconds error in RT
}
END_SECTION

START_SECTION([EXTRA_FULL_META] void store(const String& filename, MapType& map))
{
  MSExperiment exp_orig;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLSqliteHandler_1.mzML"), exp_orig);

  SqMassFile::SqMassConfig config;
  config.use_lossy_numpress = true;
  config.linear_fp_mass_acc = 0.0001;
  config.write_full_meta = true;

  {
    Precursor p;
    std::set<Precursor::ActivationMethod> tmp = {Precursor::ActivationMethod::BIRD};
    p.setActivationMethods(tmp);
    p.setActivationEnergy(500);
    p.setCharge(4);
    p.setMZ(600);
    p.setIsolationWindowUpperOffset(7);
    p.setIsolationWindowLowerOffset(14);
    p.setDriftTime(0.5);
    p.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
    p.setMetaValue("peptide_sequence", "PEPTIDEK");

    std::vector<Precursor> prec;
    prec.push_back(p);
    exp_orig.getSpectrum(0).setPrecursors(prec);
    exp_orig.getChromatogram(0).setPrecursor(p);

    Product pr;
    // pr.setCharge(6);
    pr.setMZ(300);
    pr.setIsolationWindowUpperOffset(10);
    pr.setIsolationWindowLowerOffset(15);
    std::vector<Product> prod;
    prod.push_back(pr);
    exp_orig.getSpectrum(0).setProducts(prod);
    exp_orig.getChromatogram(0).setProduct(pr);
    TEST_REAL_SIMILAR(exp_orig.getSpectrum(0).getPrecursors()[0].getActivationEnergy(), 500.0);
  }

  SqMassFile file;
  file.setConfig(config);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  std::cout << "Storing in file " << tmp_filename << std::endl;
  file.store(tmp_filename, exp_orig);

  MSExperiment exp;
  file.load(tmp_filename, exp);

  TEST_EQUAL(exp.getNrSpectra(), exp_orig.getSpectra().size())
  TEST_EQUAL(exp.getNrChromatograms(), exp_orig.getChromatograms().size())
  TEST_EQUAL(exp.getNrSpectra(), 2)
  TEST_EQUAL(exp.getNrChromatograms(), 1)
  TEST_EQUAL(exp.getSpectrum(0) == exp_orig.getSpectra()[0], false) // no exact duplicate

  {
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getRT(), exp_orig.getSpectrum(0).getRT())
    TEST_EQUAL(exp.getSpectrum(0).getNativeID(), exp_orig.getSpectrum(0).getNativeID())
    TEST_EQUAL(exp.getSpectrum(0).getMSLevel(), exp_orig.getSpectrum(0).getMSLevel())
    TEST_EQUAL(exp.getSpectrum(0).getInstrumentSettings().getPolarity(), exp_orig.getSpectrum(0).getInstrumentSettings().getPolarity())
    TEST_EQUAL(exp.getSpectrum(0).getProducts().size(), 1)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getProducts()[0].getMZ(), 300)
    // TEST_EQUAL(exp.getSpectrum(0).getProducts()[0].getCharge(), 6)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getProducts()[0].getIsolationWindowUpperOffset(), 10)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getProducts()[0].getIsolationWindowLowerOffset(), 15)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors().size(), 1)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getActivationEnergy(), 500)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors()[0].getCharge(), 4)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getDriftTime(), 0.5)
    // TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getDriftTimeUnit(), 0.5) // not implemented
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getMZ(), 600)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getIsolationWindowUpperOffset(), 7)
    TEST_REAL_SIMILAR(exp.getSpectrum(0).getPrecursors()[0].getIsolationWindowLowerOffset(), 14)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors()[0].metaValueExists("peptide_sequence"), true)
    TEST_EQUAL(exp.getSpectrum(0).getPrecursors()[0].getMetaValue("peptide_sequence"), "PEPTIDEK")
  }
  {
    TEST_EQUAL(exp.getChromatogram(0).getNativeID(), exp_orig.getChromatogram(0).getNativeID())
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getProduct().getMZ(), 300)
    // TEST_EQUAL(exp.getChromatogram(0).getProduct().getCharge(), 6)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getProduct().getIsolationWindowUpperOffset(), 10)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getProduct().getIsolationWindowLowerOffset(), 15)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getActivationEnergy(), 500)
    TEST_EQUAL(exp.getChromatogram(0).getPrecursor().getCharge(), 4)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getDriftTime(), 0.5)
    // TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getDriftTimeUnit(), 0.5) // not implemented
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getMZ(), 600)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getIsolationWindowUpperOffset(), 7)
    TEST_REAL_SIMILAR(exp.getChromatogram(0).getPrecursor().getIsolationWindowLowerOffset(), 14)
    TEST_EQUAL(exp.getChromatogram(0).getPrecursor().metaValueExists("peptide_sequence"), true)
    TEST_EQUAL(exp.getChromatogram(0).getPrecursor().getMetaValue("peptide_sequence"), "PEPTIDEK")
  }

  MSExperiment exp2 = exp_orig;

  // using full meta should give 1:1 mapping of experimental settings ...
  TEST_EQUAL(exp.getExperimentalSettings() == (OpenMS::ExperimentalSettings)exp2, true)
  TEST_EQUAL(exp.getSqlRunID(), exp2.getSqlRunID())

  // Logic of comparison: if the absolute difference criterion is fulfilled,
  // the relative one does not matter. If the absolute difference is larger
  // than allowed, the test does not fail if the relative difference is less
  // than allowed.
  // Note that the sample spectrum intensity has a very large range, from
  // 0.00013 to 183 838 intensity and encoding both values with high accuracy
  // is difficult.
  TOLERANCE_ABSOLUTE(1e-4)
  TOLERANCE_RELATIVE(1.001) // 0.1 % error for intensity 
  
  cmpDataIntensity(exp, exp2, 1e-4, 1.001);   // 0.1 % error for intensity 
  cmpDataMZ(exp, exp2, 1e-5, 1.000001);  // less than 1ppm error for m/z 
  cmpDataRT(exp, exp2, 0.05, 1.00000001);  // max 0.05 seconds error in RT
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

