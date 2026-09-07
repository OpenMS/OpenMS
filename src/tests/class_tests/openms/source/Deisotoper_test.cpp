// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>

///////////////////////////

START_TEST(Deisotoper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

		    
START_SECTION(static void deisotopeAndSingleChargeMSSpectrum(MSSpectrum& in,
                                          double fragment_tolerance, 
                                          bool fragment_unit_ppm,
                                          int min_charge = 1, 
                                          int max_charge = 3,
                                          bool keep_only_deisotoped = false,
                                          unsigned int min_isopeaks = 3, 
                                          unsigned int max_isopeaks = 10,
                                          bool make_single_charged = true,
                                          bool annotate_charge = false))
{
   MSSpectrum two_patterns;
   Peak1D p;
   p.setIntensity(1.0);

   // one charge one pattern
   p.setMZ(100.0);
   two_patterns.push_back(p);
   p.setMZ(100.0 + Constants::C13C12_MASSDIFF_U);
   two_patterns.push_back(p);
   p.setMZ(100.0 + 2.0 * Constants::C13C12_MASSDIFF_U);
   two_patterns.push_back(p);


   // one charge two pattern
   p.setMZ(200.0);
   two_patterns.push_back(p);
   p.setMZ(200.0 + 0.5 * Constants::C13C12_MASSDIFF_U);
   two_patterns.push_back(p);
   p.setMZ(200.0 + 2.0 * 0.5 * Constants::C13C12_MASSDIFF_U);
   two_patterns.push_back(p);

   MSSpectrum theo0 = two_patterns;
   Deisotoper::deisotopeAndSingleCharge(theo0, 
		   10.0, 
		   true, 
		   1, 
		   2, 
		   true, 
		   2,
		   10,
		   false, 
		   true);

   TEST_EQUAL(theo0.size(), 2); // two peaks after deisotoping
   TEST_REAL_SIMILAR(theo0[0].getMZ(), 100); 
   TEST_REAL_SIMILAR(theo0[1].getMZ(), 200); 

   theo0 = two_patterns;
   Deisotoper::deisotopeAndSingleCharge(theo0, 
		   10.0, 
		   true, 
		   1, 
		   2, 
		   true, 
		   2,
		   10,
		   true,  // convert to charge 1
		   true);

   TEST_EQUAL(theo0.size(), 2); // two peaks after deisotoping
   TEST_REAL_SIMILAR(theo0[0].getMZ(), 100); 
   TEST_REAL_SIMILAR(theo0[1].getMZ(), 400.0 - Constants::PROTON_MASS_U); 

   // create a theoretical spectrum generator 
   // and configure to add isotope patterns
   TheoreticalSpectrumGenerator spec_generator;
   Param param = spec_generator.getParameters();
   param.setValue("isotope_model", "coarse");
   param.setValue("max_isotope", 3);
   param.setValue("add_a_ions", "false");
   param.setValue("add_b_ions", "false");
   param.setValue("add_losses", "false");
   param.setValue("add_precursor_peaks", "false");
   spec_generator.setParameters(param);
   MSSpectrum theo1;
   AASequence peptide1 = AASequence::fromString("PEPTIDE");
   spec_generator.getSpectrum(theo1, peptide1, 1, 2);// charge 1..2
   TEST_EQUAL(theo1.size(), 36);
   theo1.sortByPosition();
   Deisotoper::deisotopeAndSingleCharge(theo1, 
		   10.0, 
		   true, 
		   1, 
		   2, 
		   true, 
		   2,
		   10,
		   false, 
		   true);
   // create theoretical spectrum without isotopic peaks for comparison to the deisotoped one
   param.setValue("isotope_model", "none");  // disable additional isotopes
   spec_generator.setParameters(param);
   MSSpectrum theo1_noiso;
   spec_generator.getSpectrum(theo1_noiso, peptide1, 1, 2); // charge 1..2
   TEST_EQUAL(theo1.size(), theo1_noiso.size()); // same number of peaks after deisotoping

   // load data with small intensity satellite peaks (e.g., amidation)
   MSExperiment input1;
   MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML"), input1);
   Deisotoper::deisotopeAndSingleCharge(input1[0], 
		   0.01,  // Da
		   false, 
		   1, 
		   3, 
		   false, 
		   2,
		   10,
		   false, 
		   true,
       false, // no iso peak count annotation
       true, // decreasing isotope model
       2, // enforce only starting from second peak
       true);
   std::string temp_file1 = File::getTempDirectory() + "/" + File::getUniqueName() + "_Deisotoper_output1.mzML";
   MzMLFile().store(temp_file1, input1);
   File::remove(temp_file1);

   // load data with small intensity satellite peaks (e.g., amidation)
   MSExperiment input2;
   std::string input2_path = OPENMS_GET_TEST_DATA_PATH("Deisotoper_input2.mzML");
   if (File::exists(input2_path))
   {
     MzMLFile().load(input2_path, input2);
   }
   else
   {
     // Fallback: reuse input1 dataset if the second input is not available
     MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("Deisotoper_input1.mzML"), input2);
   }
   Deisotoper::deisotopeAndSingleCharge(input2[0], 
		   0.01,  // Da
		   false, 
		   1, 
		   3, 
		   false, 
		   2,
		   10,
		   false, 
		   true,
       false, // no iso peak count annotation
       true, // decreasing isotope model
       2, // enforce only starting from second peak
       true);
   std::string temp_file2 = File::getTempDirectory() + "/" + File::getUniqueName() + "_Deisotoper_output2.mzML";
   MzMLFile().store(temp_file2, input2);
   File::remove(temp_file2);

   // Regression test for issue #10067: a precursor with unknown charge (0) must
   // not activate the precursor-mass constraint (which would use a zero mass and
   // reject every fragment cluster, emptying the spectrum with keep_only_deisotoped).
   {
     MSSpectrum base;
     Peak1D pk;
     pk.setIntensity(1.0);
     // one doubly-charged isotope cluster (0.5 Da spacing)
     pk.setMZ(200.0);                                          base.push_back(pk);
     pk.setMZ(200.0 + 0.5 * Constants::C13C12_MASSDIFF_U);     base.push_back(pk);
     pk.setMZ(200.0 + 1.0 * Constants::C13C12_MASSDIFF_U);     base.push_back(pk);

     // (a) no precursor: cluster is detected and collapsed to its monoisotopic peak
     MSSpectrum s = base;
     Deisotoper::deisotopeAndSingleCharge(s, 10.0, true, 1, 2, true, 2, 10, false, true);
     TEST_EQUAL(s.size(), 1);
     const MSSpectrum without_precursor = s;

     // (b) precursor present but charge unknown (0): must behave like (a), not empty the spectrum
     s = base;
     Precursor prec_unknown;
     prec_unknown.setMZ(200.0);
     prec_unknown.setCharge(0);
     s.setPrecursors({prec_unknown});
     Deisotoper::deisotopeAndSingleCharge(s, 10.0, true, 1, 2, true, 2, 10, false, true);
     TEST_EQUAL(s.size(), without_precursor.size());
     TEST_REAL_SIMILAR(s[0].getMZ(), without_precursor[0].getMZ());
     TEST_REAL_SIMILAR(s[0].getMZ(), 200.0);
     TEST_TRUE(s.getIntegerDataArrays() == without_precursor.getIntegerDataArrays());
     TEST_EQUAL(s.getIntegerDataArrays()[0][0], 2);

     // (c) known precursor charge large enough to keep the fragment cluster: cluster retained
     //     and annotated with the detected charge (2)
     s = base;
     Precursor prec_known;
     prec_known.setMZ(2000.0);
     prec_known.setCharge(2);
     s.setPrecursors({prec_known});
     Deisotoper::deisotopeAndSingleCharge(s, 10.0, true, 1, 2, true, 2, 10, false, true);
     TEST_EQUAL(s.size(), 1);
     TEST_EQUAL(s.getIntegerDataArrays().size(), 1);
     TEST_EQUAL(s.getIntegerDataArrays()[0].size(), 1);
     TEST_EQUAL(s.getIntegerDataArrays()[0][0], 2);

     // (d) known precursor charge whose neutral mass is below the fragment cluster: the
     //     constraint is still functional and rejects the cluster (empty with keep_only)
     s = base;
     Precursor prec_low;
     prec_low.setMZ(100.0);
     prec_low.setCharge(1);
     s.setPrecursors({prec_low});
     Deisotoper::deisotopeAndSingleCharge(s, 10.0, true, 1, 2, true, 2, 10, false, true);
     TEST_EQUAL(s.size(), 0);

     // (e) retention mode (keep_only_deisotoped=false) with an unassigned peak and unknown
     //     precursor charge: the cluster collapses to its monoisotopic peak and the
     //     unrelated peak is retained rather than dropped
     MSSpectrum base_plus = base;
     Peak1D lone;
     lone.setIntensity(1.0);
     lone.setMZ(600.0);
     base_plus.push_back(lone);
     base_plus.sortByPosition();
     MSSpectrum retained_without_precursor = base_plus;
     Deisotoper::deisotopeAndSingleCharge(retained_without_precursor, 10.0, true, 1, 2, false, 2, 10, false, true);
     s = base_plus;
     s.setPrecursors({prec_unknown});
     Deisotoper::deisotopeAndSingleCharge(s, 10.0, true, 1, 2, false, 2, 10, false, true);
     TEST_EQUAL(s.size(), 2);
     TEST_REAL_SIMILAR(s[0].getMZ(), 200.0);
     TEST_REAL_SIMILAR(s[1].getMZ(), 600.0);
     TEST_EQUAL(s.size(), retained_without_precursor.size());
     TEST_REAL_SIMILAR(s[0].getMZ(), retained_without_precursor[0].getMZ());
     TEST_REAL_SIMILAR(s[1].getMZ(), retained_without_precursor[1].getMZ());
     TEST_TRUE(s.getIntegerDataArrays() == retained_without_precursor.getIntegerDataArrays());
     TEST_EQUAL(s.getIntegerDataArrays()[0][0], 2);
     TEST_EQUAL(s.getIntegerDataArrays()[0][1], 0);
   }
}
END_SECTION

START_SECTION(static void deisotopeWithAveragineModel(MSSpectrum& spectrum,
                                                      double fragment_tolerance,
                                                      bool fragment_unit_ppm))
{
  // spectrum with one isotopic pattern
  MSSpectrum spec;
  CoarseIsotopePatternGenerator gen(5);
  IsotopeDistribution distr = gen.estimateFromPeptideWeight(700);
  double base_mz1 = distr[0].getMZ();
  for (auto it = distr.begin(); it != distr.end(); ++it)
  {
    if (it->getIntensity() != 0)
    {
      it->setIntensity(it->getIntensity() * 10);
      spec.push_back(*it);
    }
  }
  spec.sortByPosition();
  MSSpectrum theo(spec);
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true);
  TEST_EQUAL(theo.size(), 1);
  TEST_REAL_SIMILAR(theo[0].getMZ(), base_mz1);

  // Peaks before and after spectrum should not be chosen
  // Shows a fault of the old algorithm occurring with e.g. deamination
  double correct_monoiso = theo.getBasePeak()->getMZ();
  double deamin_mz = spec.front().getMZ() - OpenMS::Constants::NEUTRON_MASS_U;
  Peak1D deamin_peak = Peak1D(deamin_mz, 0.06f);
  spec.push_back(deamin_peak);
  spec.sortByPosition();
  theo = spec;
  MSSpectrum theo1(spec);
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true, 5000, 1, 3, true);//keep only deisotoped
  TEST_REAL_SIMILAR(theo.front().getMZ(), correct_monoiso);
  Deisotoper::deisotopeAndSingleCharge(theo1, 10.0, true, 1, 3, true);
  TEST_NOT_EQUAL(theo1.front().getMZ(), correct_monoiso);// passes -> not equal

  // Test a peak with zero intensitiy
  double add_mz = spec.back().getMZ() + OpenMS::Constants::C13C12_MASSDIFF_U;
  Peak1D add_peak(add_mz, 0.0);
  spec.push_back(add_peak);
  theo = spec;
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true);
  TEST_NOT_EQUAL(theo.back().getIntensity(), 0.0);// the new peak should be removed

  // Additional peaks that only fit m/z - wise should not disturb cluster formation
  spec.back().setIntensity(20);// intensity is a lot too high to fit correct distribution
  theo = spec;
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true, -1);// do not remove low intensities
  TEST_EQUAL(theo.size(), 3);
  TEST_REAL_SIMILAR(theo.back().getMZ(), add_mz);// last peak is still there

  // spectrum with two isotopic patterns
  distr = gen.estimateFromPeptideWeight(500);
  double base_mz2 = distr[0].getMZ();
  for (auto it = distr.begin(); it != distr.end(); ++it)
  {
    if (it->getIntensity() != 0)
    {
      it->setMZ((it->getMZ() + OpenMS::Constants::PROTON_MASS_U) / 2);// set to charge 2
      spec.push_back((*it));
    }
  }
  theo = spec;
  theo.sortByPosition();
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true, 5000, 1, 3, true);// keep only deisotoped
  TEST_EQUAL(theo.size(), 2);
  TEST_EQUAL(theo[0].getMZ(), base_mz2);
  TEST_EQUAL(theo[1].getMZ(), base_mz1);

  // Add unassignable peaks
  Peak1D peak1(550, 0.8f);
  spec.push_back(peak1);
  Peak1D peak2(600, 0.9f);
  spec.push_back(peak2);
  spec.sortByPosition();
  theo = spec;
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true, -1);// do not remove low intensities
  TEST_EQUAL(theo.size(), 6);                                      // two spectra, one peak before, one after one spectrum, and two unassignable peaks

  // keep only deisotoped
  theo = spec;
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true, 5000, 1, 3, true); // keep only deisotoped
  TEST_EQUAL(theo.size(), 2);

  // test with complete theoretical spectrum

  // create a theoretical spectrum generator
  // and configure to add isotope patterns
  TheoreticalSpectrumGenerator spec_generator;
  Param param = spec_generator.getParameters();
  param.setValue("isotope_model", "coarse");
  param.setValue("max_isotope", 3);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "false");
  param.setValue("add_losses", "false");
  param.setValue("add_precursor_peaks", "false");
  spec_generator.setParameters(param);
  param.setValue("isotope_model", "coarse");
  spec_generator.setParameters(param);

  AASequence peptide1 = AASequence::fromString("PEPTIDE");

  theo.clear(true);
  spec_generator.getSpectrum(theo, peptide1, 1, 2);// charge 1..2
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true);

  // create theoretical spectrum without isotopic peaks for comparison to the deisotoped one
  param.setValue("isotope_model", "none");// disable additional isotopes
  spec_generator.setParameters(param);
  MSSpectrum theo_noiso;
  spec_generator.getSpectrum(theo_noiso, peptide1, 1, 2);// charge 1..2
  TEST_EQUAL(theo.size(), theo_noiso.size());            // same number of peaks after deisotoping

  // simpler tests with patterns where all isotopic peaks have the same intensity
  MSSpectrum two_patterns;
  Peak1D p;
  two_patterns.clear(true);
  p.setIntensity(1.0);

  // first pattern
  p.setMZ(100.0);
  two_patterns.push_back(p);
  p.setMZ(100.0 + Constants::C13C12_MASSDIFF_U);
  two_patterns.push_back(p);
  p.setMZ(100.0 + 2.0 * Constants::C13C12_MASSDIFF_U);
  two_patterns.push_back(p);

  // second pattern
  p.setMZ(200.0);
  two_patterns.push_back(p);
  p.setMZ(200.0 + 0.5 * Constants::C13C12_MASSDIFF_U);
  two_patterns.push_back(p);
  p.setMZ(200.0 + 2.0 * 0.5 * Constants::C13C12_MASSDIFF_U);
  two_patterns.push_back(p);
  theo = two_patterns;
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true);
  TEST_EQUAL(theo.size(), 6);// all six peaks remain, since the patterns should not be similar to averagine model

  // Test with a section of an actual spectrum
  MzMLFile file;
  PeakMap exp;
  file.load(OPENMS_GET_TEST_DATA_PATH("Deisotoper_test_in.mzML"), exp);
  theo.clear(true);
  theo = exp.getSpectrum(0);// copy for readability
  theo1.clear(true);
  theo1 = exp.getSpectrum(0);// for next test
  Size ori_size = theo.size();
  Deisotoper::deisotopeWithAveragineModel(theo, 10.0, true, 5000, 1, 3, true);// keep only deisotoped
  TEST_NOT_EQUAL(theo.size(), ori_size);
  file.load(OPENMS_GET_TEST_DATA_PATH("Deisotoper_test_out.mzML"), exp);
  TEST_EQUAL(theo, exp.getSpectrum(0));

  // Test if the algorithm also works if we do not remove the low (and zero) intensity peaks
  Deisotoper::deisotopeWithAveragineModel(theo1, 10.0, true, -1, 1, 3, true);// do not remove low intensity peaks beforehand, but keep only deisotoped
  TEST_EQUAL(theo1.size(), 104);

  // Regression test for issue #10067: a precursor with unknown charge (0) must
  // not activate the precursor-mass constraint (which would use a zero mass and
  // reject every fragment cluster, emptying the spectrum with keep_only_deisotoped).
  {
    // build a single doubly-charged isotope cluster (~m/z 500) that matches the averagine model
    CoarseIsotopePatternGenerator cluster_gen(5);
    IsotopeDistribution cluster_distr = cluster_gen.estimateFromPeptideWeight(1000);
    MSSpectrum base;
    for (auto it = cluster_distr.begin(); it != cluster_distr.end(); ++it)
    {
      if (it->getIntensity() != 0)
      {
        Peak1D pk;
        pk.setMZ((it->getMZ() + Constants::PROTON_MASS_U) / 2.0);// charge 2
        pk.setIntensity(it->getIntensity() * 10);
        base.push_back(pk);
      }
    }
    base.sortByPosition();

    // (a) no precursor: cluster is detected and collapsed to its monoisotopic peak
    MSSpectrum s = base;
    Deisotoper::deisotopeWithAveragineModel(s, 10.0, true, 5000, 1, 3, true, 2, 10, true, true);// annotate charge
    TEST_EQUAL(s.size(), 1);
    const MSSpectrum without_precursor = s;

    // (b) precursor present but charge unknown (0): must behave like (a), not empty the spectrum
    s = base;
    Precursor prec_unknown;
    prec_unknown.setMZ(500.0);
    prec_unknown.setCharge(0);
    s.setPrecursors({prec_unknown});
    Deisotoper::deisotopeWithAveragineModel(s, 10.0, true, 5000, 1, 3, true, 2, 10, true, true);
    TEST_EQUAL(s.size(), without_precursor.size());
    TEST_REAL_SIMILAR(s[0].getMZ(), without_precursor[0].getMZ());
    TEST_TRUE(s.getIntegerDataArrays() == without_precursor.getIntegerDataArrays());
    TEST_EQUAL(s.getIntegerDataArrays()[0][0], 2);

    // (c) known precursor charge large enough to keep the fragment cluster: cluster retained
    //     and annotated with the detected charge (2)
    s = base;
    Precursor prec_known;
    prec_known.setMZ(2000.0);
    prec_known.setCharge(2);
    s.setPrecursors({prec_known});
    Deisotoper::deisotopeWithAveragineModel(s, 10.0, true, 5000, 1, 3, true, 2, 10, true, true);// annotate charge
    TEST_EQUAL(s.size(), 1);
    TEST_EQUAL(s.getIntegerDataArrays().size(), 1);
    TEST_EQUAL(s.getIntegerDataArrays()[0].size(), 1);
    TEST_EQUAL(s.getIntegerDataArrays()[0][0], 2);

    // (d) known precursor charge whose neutral mass is below the fragment cluster: the
    //     constraint is still functional and rejects the cluster (empty with keep_only)
    s = base;
    Precursor prec_low;
    prec_low.setMZ(100.0);
    prec_low.setCharge(1);
    s.setPrecursors({prec_low});
    Deisotoper::deisotopeWithAveragineModel(s, 10.0, true, 5000, 1, 3, true);
    TEST_EQUAL(s.size(), 0);

    // (e) retention mode (keep_only_deisotoped=false) with an unassigned peak and unknown
    //     precursor charge: the unrelated peak is retained rather than dropped
    MSSpectrum base_plus = base;
    Peak1D lone;
    lone.setIntensity(50.0);
    lone.setMZ(800.0);
    base_plus.push_back(lone);
    base_plus.sortByPosition();
    MSSpectrum retained_without_precursor = base_plus;
    Deisotoper::deisotopeWithAveragineModel(retained_without_precursor, 10.0, true, -1, 1, 3, false, 2, 10, true, true);
    s = base_plus;
    s.setPrecursors({prec_unknown});
    Deisotoper::deisotopeWithAveragineModel(s, 10.0, true, -1, 1, 3, false, 2, 10, true, true);// keep unassigned peaks and annotate charge
    TEST_EQUAL(s.size(), 2);
    TEST_REAL_SIMILAR(s[0].getMZ(), 800.0);// the unassigned peak survived (sorts before the converted monoisotopic peak)
    TEST_EQUAL(s.size(), retained_without_precursor.size());
    TEST_REAL_SIMILAR(s[0].getMZ(), retained_without_precursor[0].getMZ());
    TEST_REAL_SIMILAR(s[1].getMZ(), retained_without_precursor[1].getMZ());
    TEST_TRUE(s.getIntegerDataArrays() == retained_without_precursor.getIntegerDataArrays());
    TEST_EQUAL(s.getIntegerDataArrays()[0][0], 0);
    TEST_EQUAL(s.getIntegerDataArrays()[0][1], 2);
  }
}
END_SECTION

START_SECTION((static bool isToleranceSupported(double fragment_tolerance, bool fragment_unit_ppm)))
{
  // Non-throwing mirror of the deisotopeAndSingleCharge() precondition (<= 100 ppm / <= 0.1 Da).
  // ppm: supported iff tolerance <= 100
  TEST_EQUAL(Deisotoper::isToleranceSupported(100.0, true), true)
  TEST_EQUAL(Deisotoper::isToleranceSupported(20.0, true), true)
  TEST_EQUAL(Deisotoper::isToleranceSupported(100.0001, true), false)
  TEST_EQUAL(Deisotoper::isToleranceSupported(150.0, true), false)
  // Da: supported iff tolerance <= 0.1
  TEST_EQUAL(Deisotoper::isToleranceSupported(0.1, false), true)
  TEST_EQUAL(Deisotoper::isToleranceSupported(0.02, false), true)
  TEST_EQUAL(Deisotoper::isToleranceSupported(0.10001, false), false)
  TEST_EQUAL(Deisotoper::isToleranceSupported(0.5, false), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION((regression: precursor mass constraint uses atomic mass units for unequal charges))
{
  // A precursor at m/z 500 with charge 2 has neutral mass 997.985447 Da.
  // Singly charged fragments at m/z 998.5 and 999.5 lie below and above
  // that mass, respectively. Using the proton mass in kg accepts both.
  Precursor precursor;
  precursor.setMZ(500.0);
  precursor.setCharge(2);

  for (const double fragment_mz : {998.5, 999.5})
  {
    MSSpectrum base;
    for (Size i = 0; i < 3; ++i)
    {
      Peak1D peak;
      peak.setMZ(fragment_mz + i * Constants::C13C12_MASSDIFF_U);
      peak.setIntensity(100.0 / (i + 1));
      base.push_back(peak);
    }

    // Ensure the isotope pattern itself is detected without a precursor limit.
    MSSpectrum unconstrained = base;
    Deisotoper::deisotopeAndSingleCharge(unconstrained, 10.0, true, 1, 1, true, 3, 10, false, true);
    TEST_EQUAL(unconstrained.size(), 1);
    TEST_REAL_SIMILAR(unconstrained[0].getMZ(), fragment_mz);
    TEST_EQUAL(unconstrained.getIntegerDataArrays()[0][0], 1);

    for (const bool fragment_unit_ppm : {false, true})
    {
      for (const bool keep_only_deisotoped : {false, true})
      {
        MSSpectrum spectrum = base;
        spectrum.setPrecursors({precursor});
        Deisotoper::deisotopeAndSingleCharge(spectrum,
          fragment_unit_ppm ? 10.0 : 0.01, fragment_unit_ppm,
          1, 1, keep_only_deisotoped, 3, 10, false, true);

        const bool below_precursor_mass = fragment_mz == 998.5;
        const Size expected_size = below_precursor_mass ? 1 : (keep_only_deisotoped ? 0 : 3);
        TEST_EQUAL(spectrum.size(), expected_size);
        TEST_EQUAL(spectrum.getIntegerDataArrays().size(), 1);
        TEST_EQUAL(spectrum.getIntegerDataArrays()[0].getName(), "charge");
        TEST_EQUAL(spectrum.getIntegerDataArrays()[0].size(), expected_size);
        for (Size i = 0; i < spectrum.size(); ++i)
        {
          TEST_REAL_SIMILAR(spectrum[i].getMZ(), base[i].getMZ());
          TEST_EQUAL(spectrum.getIntegerDataArrays()[0][i], below_precursor_mass ? 1 : 0);
        }
      }
    }
  }
}
END_SECTION

END_TEST
