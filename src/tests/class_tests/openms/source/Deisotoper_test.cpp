// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <iostream>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

START_TEST(Deisotoper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


START_SECTION(static void deisotopeAndSingleChargeMSSpectrum(MSSpectrum& in,
                                          double fragment_tolerance,
                                          bool fragment_unit_ppm,
                                          std::string model = "none",
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
       "decreasing",
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
       "decreasing",
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
   param.setValue("add_isotopes", "true");
   param.setValue("max_isotope", 3);
   param.setValue("add_a_ions", "false");
   param.setValue("add_b_ions", "false");
   param.setValue("add_losses", "false");
   param.setValue("add_precursor_peaks", "false");
   spec_generator.setParameters(param);
   MSSpectrum theo1;
   AASequence peptide1 = AASequence::fromString("PEPTIDE");
   spec_generator.getSpectrum(theo1, peptide1, 1, 2); // charge 1..2
   TEST_EQUAL(theo1.size(), 36);
   theo1.sortByPosition();
   Deisotoper::deisotopeAndSingleCharge(theo1,
       10.0,
       true,
       "decreasing",
       1,
       2,
       true,
       2,
       10,
       false,
       true);
   // create theoretical spectrum without isotopic peaks for comparison to the deisotoped one
   param.setValue("add_isotopes", "false");  // disable additional isotopes
   spec_generator.setParameters(param);
   MSSpectrum theo1_noiso;
   spec_generator.getSpectrum(theo1_noiso, peptide1, 1, 2); // charge 1..2
   TEST_EQUAL(theo1.size(), theo1_noiso.size()); // same number of peaks after deisotoping

   // test peptide DFPIANGER (mono-isotopic peak is highest)
   {
     int z = 5;
     auto seq = AASequence::fromString("DFPIANGER");
     auto seq_formula = seq.getFormula() + EmpiricalFormula("H" + String(z));
     auto isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(5) );

     MSSpectrum s;
     for (auto iso : isotopes.getContainer() )
     {
       iso.setMZ( iso.getMZ() / z );
       s.push_back(iso);
     }
     Deisotoper::deisotopeAndSingleCharge(s,
         10.0, true,
         "decreasing",
         1, z,
         true,
         2, 10,
         true,
         true);
     TEST_EQUAL(s.size(), 1);
     // expect a single peak at the location of the [M+H]+ ion
     TEST_REAL_SIMILAR(s[0].getMZ(), seq.getMonoWeight(Residue::ResidueType::Full, 1));
     TEST_REAL_SIMILAR(s[0].getMZ(), 1018.49798342972);
     // test charge estimation
     TEST_EQUAL(s.getIntegerDataArrays()[0][0], z);
   }

   // test peptide DFPIANGERDFPIANGERDFPIANGERDFPIANGER (mono-isotopic peak is not the highest)
   {
     int z = 5;
     auto seq = AASequence::fromString("DFPIANGERDFPIANGERDFPIANGERDFPIANGER");
     auto seq_formula = seq.getFormula() + EmpiricalFormula("H" + String(z));
     auto isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(8) );

     MSSpectrum s;
     for (auto iso : isotopes.getContainer() ) // copy data
     {
       iso.setMZ( iso.getMZ() / z );
       s.push_back(iso);
     }
     Deisotoper::deisotopeAndSingleCharge(s,
         10.0, true,
         "none",
         1, z,
         true,
         2, 10,
         true,
         true);
     // expect a single peak at the location of the [M+H]+ ion
     TEST_EQUAL(s.size(), 1);
     TEST_REAL_SIMILAR(s[0].getMZ(), seq.getMonoWeight(Residue::ResidueType::Full, 1));
     // test charge estimation
     TEST_EQUAL(s.getIntegerDataArrays()[0][0], z);
   }

   // test peptide DFPIANGERDFPIANGERDFPIANGERDFPIANGER (mono-isotopic peak is not the highest) with wrong settings
   {
     int z = 5;
     auto seq = AASequence::fromString("DFPIANGERDFPIANGERDFPIANGERDFPIANGER");
     auto seq_formula = seq.getFormula() + EmpiricalFormula("H" + String(z));
     IsotopeDistribution isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(8) );

     MSSpectrum s;
     for (auto iso : isotopes.getContainer() )
     {
       iso.setMZ( iso.getMZ() / z );
       s.push_back(iso);
     }
     Deisotoper::deisotopeAndSingleCharge(s,
         10.0, true,
         "decreasing", // WRONG - wont work
         1, z,
         true,
         2, 10,
         true,
         true);

     // WRONG: model creates 3 peaks instead of a single peak
     TEST_EQUAL(s.size(), 3);
     TEST_REAL_SIMILAR(s[0].getMZ(), isotopes.getContainer()[0].getMZ()/z ); // wrong peak / left in place
     TEST_REAL_SIMILAR(s[1].getMZ(), isotopes.getContainer()[1].getMZ()/z ); // wrong peak / left in place
     TEST_REAL_SIMILAR(s[2].getMZ(), seq.getMonoWeight(Residue::ResidueType::Full, 3)); // 3rd isotope!
     // test charge estimation
     TEST_EQUAL(s.getIntegerDataArrays()[0][2], z);
   }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
