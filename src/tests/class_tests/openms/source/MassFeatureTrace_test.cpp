// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassFeatureTrace, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassFeatureTrace* ptr = 0;
MassFeatureTrace* null_ptr = 0;
START_SECTION(MassFeatureTrace())
{
  ptr = new MassFeatureTrace();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MassFeatureTrace())
{
  delete ptr;
}
END_SECTION

/// sample input for testing ///

MassFeatureTrace mass_tracer;

MSSpectrum sample_spec;
sample_spec.setRT(50.0);
sample_spec.setMSLevel(1);
DeconvolvedSpectrum deconv_spec1(1);
deconv_spec1.setOriginalSpectrum(sample_spec);

PeakGroup tmp_pg = PeakGroup(15, 18, true);
auto p1 = new Peak1D(1000.8455675085044, 8347717.5);
FLASHDeconvHelperStructs::LogMzPeak tmp_p1(*p1, true);
tmp_p1.abs_charge = 18;
tmp_p1.isotopeIndex = 8;

p1 = new Peak1D(1000.9013094439375, 10087364);
FLASHDeconvHelperStructs::LogMzPeak tmp_p2(*p1, true);
tmp_p2.abs_charge = 18;
tmp_p2.isotopeIndex = 9;

p1 = new Peak1D(1000.9570513793709, 11094268);
FLASHDeconvHelperStructs::LogMzPeak tmp_p3(*p1, true);
tmp_p3.abs_charge = 18;
tmp_p3.isotopeIndex = 10;

p1 = new Peak1D(1001.0127933148044, 11212854);
FLASHDeconvHelperStructs::LogMzPeak tmp_p4(*p1, true);
tmp_p4.abs_charge = 18;
tmp_p4.isotopeIndex = 11;

p1 = new Peak1D(1001.0685352502376, 10497022);
FLASHDeconvHelperStructs::LogMzPeak tmp_p5(*p1, true);
tmp_p5.abs_charge = 18;
tmp_p5.isotopeIndex = 12;

p1 = new Peak1D(1001.124277185671, 9162559);
FLASHDeconvHelperStructs::LogMzPeak tmp_p6(*p1, true);
tmp_p6.abs_charge = 18;
tmp_p6.isotopeIndex = 13;

p1 = new Peak1D(1059.6595846286061, 8347717.5);
FLASHDeconvHelperStructs::LogMzPeak tmp_p7(*p1, true);
tmp_p7.abs_charge = 17;
tmp_p7.isotopeIndex = 8;

p1 = new Peak1D(1059.7186055014179, 10087364);
FLASHDeconvHelperStructs::LogMzPeak tmp_p8(*p1, true);
tmp_p8.abs_charge = 17;
tmp_p8.isotopeIndex = 9;

p1 = new Peak1D(1059.7776263742296, 11094268);
FLASHDeconvHelperStructs::LogMzPeak tmp_p9(*p1, true);
tmp_p9.abs_charge = 17;
tmp_p9.isotopeIndex = 10;

p1 = new Peak1D(1059.8366472470416, 11212854);
FLASHDeconvHelperStructs::LogMzPeak tmp_p10(*p1, true);
tmp_p10.abs_charge = 17;
tmp_p10.isotopeIndex = 11;

p1 = new Peak1D(1059.8956681198531, 10497022);
FLASHDeconvHelperStructs::LogMzPeak tmp_p11(*p1, true);
tmp_p11.abs_charge = 17;
tmp_p11.isotopeIndex = 12;

p1 = new Peak1D(1059.9546889926651, 9162559);
FLASHDeconvHelperStructs::LogMzPeak tmp_p12(*p1, true);
tmp_p12.abs_charge = 17;
tmp_p12.isotopeIndex = 13;

tmp_pg.push_back(tmp_p1);
tmp_pg.push_back(tmp_p2);
tmp_pg.push_back(tmp_p3);
tmp_pg.push_back(tmp_p4);
tmp_pg.push_back(tmp_p5);
tmp_pg.push_back(tmp_p6);
tmp_pg.push_back(tmp_p7);
tmp_pg.push_back(tmp_p8);
tmp_pg.push_back(tmp_p9);
tmp_pg.push_back(tmp_p10);
tmp_pg.push_back(tmp_p11);
tmp_pg.push_back(tmp_p12);
tmp_pg.updateMonomassAndIsotopeIntensities();
deconv_spec1.push_back(tmp_pg);

sample_spec.setRT(55.0);
DeconvolvedSpectrum deconv_spec2(2);
deconv_spec2.setOriginalSpectrum(sample_spec);
deconv_spec2.push_back(tmp_pg);

sample_spec.setRT(61.0);
DeconvolvedSpectrum deconv_spec3(3);
deconv_spec3.setOriginalSpectrum(sample_spec);
deconv_spec3.push_back(tmp_pg);
//////////////////////////////

mass_tracer.storeInformationFromDeconvolvedSpectrum(deconv_spec1);
mass_tracer.storeInformationFromDeconvolvedSpectrum(deconv_spec2);
mass_tracer.storeInformationFromDeconvolvedSpectrum(deconv_spec3);

/// < public methods without tests >
/// - storeInformationFromDeconvolvedSpectrum : only private variables are affected (cannot test)
/// - copy, assignment, move constructor -> not used.
/// - size
START_SECTION((std::vector<FLASHDeconvHelperStructs::MassFeature> findFeatures(const PrecalculatedAveragine &averagine)))
{
  // prepare findFeature arguments
  std::unordered_map<int, PeakGroup> null_map;
  FLASHDeconvAlgorithm fd = FLASHDeconvAlgorithm();
  Param fd_param;
  fd_param.setValue("min_charge", 5);
  fd_param.setValue("max_charge", 20);
  fd_param.setValue("max_mass", 50000.);
  fd.setParameters(fd_param);
  fd.calculateAveragine(false);
  FLASHDeconvHelperStructs::PrecalculatedAveragine averagine = fd.getAveragine();

  std::vector<FLASHDeconvHelperStructs::MassFeature> found_feature = mass_tracer.findFeatures(averagine);
  OPENMS_LOG_INFO << found_feature.size() << std::endl;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


