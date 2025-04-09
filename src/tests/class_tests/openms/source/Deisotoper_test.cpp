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
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
