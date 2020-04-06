// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <iostream>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>

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
   spec_generator.getSpectrum(theo1, peptide1, 1, 2); // charge 1..2
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
