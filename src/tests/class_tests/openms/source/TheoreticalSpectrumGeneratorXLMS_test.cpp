// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <iostream>


START_TEST(TheoreticalSpectrumGeneratorXLMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TheoreticalSpectrumGeneratorXLMS* ptr = nullptr;
TheoreticalSpectrumGeneratorXLMS* nullPointer = nullptr;

/// mostly copied from TheoreticalSpectrumGenerator_test
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
START_SECTION(TheoreticalSpectrumGeneratorXLMS())
  ptr = new TheoreticalSpectrumGeneratorXLMS();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(TheoreticalSpectrumGeneratorXLMS(const TheoreticalSpectrumGeneratorXLMS& source))
  TheoreticalSpectrumGeneratorXLMS copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~TheoreticalSpectrumGeneratorXLMS())
  delete ptr;
END_SECTION

ptr = new TheoreticalSpectrumGeneratorXLMS();
AASequence peptide = AASequence::fromString("IFSQVGK");

START_SECTION(TheoreticalSpectrumGeneratorXLMS& operator = (const TheoreticalSpectrumGeneratorXLMS& tsg))
  TheoreticalSpectrumGeneratorXLMS copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

START_SECTION(virtual void getCommonIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, bool frag_alpha, int charge = 1, Size link_pos_2 = 0))
  PeakSpectrum spec;
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 2);
  TEST_EQUAL(spec.size(), 12)

  TOLERANCE_ABSOLUTE(0.001)

  double result[] = {57.54930, 74.06004, 102.57077, 114.09134, 131.08351, 147.11280, 152.10497, 174.59953, 204.13426, 261.15975, 303.20268, 348.19178};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 3);
  TEST_EQUAL(spec.size(), 18)

  spec.clear(true);
  Param param(ptr->getParameters());
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 3);
  TEST_EQUAL(spec.size(), 54)


//  // test annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 3);

  // 6 ion types with 3 charges each are expected
  TEST_EQUAL(spec.size(), 18)

  set<String> ion_names;
  ion_names.insert("[alpha|ci$b1]");
  ion_names.insert("[alpha|ci$b2]");
  ion_names.insert("[alpha|ci$b3]");
  ion_names.insert("[alpha|ci$x1]");
  ion_names.insert("[alpha|ci$x2]");
  ion_names.insert("[alpha|ci$x3]");

  PeakSpectrum::StringDataArray string_array = spec.getStringDataArrays().at(0);

  // check if all ion names have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    String name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // beta annotations
  spec.clear(true);
  ptr->getCommonIonSpectrum(spec, peptide, 3, false, 3);
  ion_names.clear();
  ion_names.insert("[beta|ci$b1]");
  ion_names.insert("[beta|ci$b2]");
  ion_names.insert("[beta|ci$b3]");
  ion_names.insert("[beta|ci$x1]");
  ion_names.insert("[beta|ci$x2]");
  ion_names.insert("[beta|ci$x3]");

  string_array = spec.getStringDataArrays().at(0);

  for (Size i = 0; i != spec.size(); ++i)
  {
    String name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // test for charges stored in IntegerDataArray
  PeakSpectrum::IntegerDataArray charge_array = spec.getIntegerDataArrays().at(0);

  int charge_counts[3] = {0, 0, 0};
  for (Size i = 0; i != spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 6)
  TEST_EQUAL(charge_counts[1], 6)
  TEST_EQUAL(charge_counts[2], 6)

  // the smallest examples, that make sense for cross-linking
  spec.clear(true);
  ptr->getCommonIonSpectrum(spec, AASequence::fromString("HA"), 0, true, 1);
  TEST_EQUAL(spec.size(), 1)

  spec.clear(true);
  ptr->getCommonIonSpectrum(spec, AASequence::fromString("HA"), 1, true, 1);
  TEST_EQUAL(spec.size(), 1)

  // loop link
  spec.clear(true);
  ptr->getCommonIonSpectrum(spec, AASequence::fromString("PEPTIDESAREWEIRD"), 1, true, 1, 14);
  TEST_EQUAL(spec.size(), 2)

  spec.clear(true);
  ptr->getCommonIonSpectrum(spec, AASequence::fromString("PEPTIDESAREWEIRD"), 2, false, 1, 14);
  TEST_EQUAL(spec.size(), 3)

  // test isotopic peaks
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 1); // not very useful combination, but it should at least run
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 3);
  // 6 ion types with 3 charges each are expected
  TEST_EQUAL(spec.size(), 18)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 2); //
  ptr->setParameters(param);
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 3);
  // 6 ion types with 3 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 36)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 3); // not supported yet, but it should at least run (with the maximal possible number of peaks)
  ptr->setParameters(param);
  ptr->getCommonIonSpectrum(spec, peptide, 3, true, 3);
  // 6 ion types with 3 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 36)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  // for quick benchmarking of implementation chances
//  param = ptr->getParameters();
//  param.setValue("add_a_ions", "true");
//  param.setValue("add_b_ions", "true");
//  param.setValue("add_c_ions", "true");
//  param.setValue("add_x_ions", "true");
//  param.setValue("add_y_ions", "true");
//  param.setValue("add_z_ions", "true");
//  param.setValue("add_metainfo", "true");
//  ptr->setParameters(param);
//  AASequence tmp_peptide = AASequence::fromString("PEPTIDEPEPTIDEPEPTIDE");
//  for (Size i = 0; i != 1e4; ++i)
//  {
//    PeakSpectrum spec;
//    ptr->getCommonIonSpectrum(spec, tmp_peptide, 9, true, 5);
//  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

END_SECTION

START_SECTION(virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2 = 0))

  // reinitialize TSG to standard parameters
  Param param(ptr->getParameters());
  param.setValue("add_isotopes", "false");
  param.setValue("max_isotope", 2);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);

  PeakSpectrum spec;
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 3);
  TEST_EQUAL(spec.size(), 12)

  TOLERANCE_ABSOLUTE(0.001)

  double result[] = {551.94577, 566.94214, 580.95645, 599.96494, 618.97210, 629.97925, 827.41502, 849.90957, 870.93103, 899.44378, 927.95451, 944.46524};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 4);
  TEST_EQUAL(spec.size(), 18)

  spec.clear(true);
//  Param param(ptr->getParameters());
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 4);
  TEST_EQUAL(spec.size(), 54)


//  // test annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "true");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);

  // 6 ion types with 4 charges each are expected
  TEST_EQUAL(spec.size(), 24)

  set<String> ion_names;
  ion_names.insert("[alpha|xi$b4]");
  ion_names.insert("[alpha|xi$b5]");
  ion_names.insert("[alpha|xi$b6]");
  ion_names.insert("[alpha|xi$x4]");
  ion_names.insert("[alpha|xi$x5]");
  ion_names.insert("[alpha|xi$x6]");

  PeakSpectrum::StringDataArray string_array = spec.getStringDataArrays().at(0);

  // check if all ion names have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    String name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // beta annotations
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, false, 2, 4);
  ion_names.clear();
  ion_names.insert("[beta|xi$b4]");
  ion_names.insert("[beta|xi$b5]");
  ion_names.insert("[beta|xi$b6]");
  ion_names.insert("[beta|xi$x4]");
  ion_names.insert("[beta|xi$x5]");
  ion_names.insert("[beta|xi$x6]");

  string_array = spec.getStringDataArrays().at(0);

  for (Size i = 0; i != spec.size(); ++i)
  {
    String name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // test for charges stored in IntegerDataArray
  PeakSpectrum::IntegerDataArray charge_array = spec.getIntegerDataArrays().at(0);

  int charge_counts[5] = {0, 0, 0, 0, 0};
  for (Size i = 0; i != spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 0)
  TEST_EQUAL(charge_counts[1], 6)
  TEST_EQUAL(charge_counts[2], 6)
  TEST_EQUAL(charge_counts[3], 6)
  TEST_EQUAL(charge_counts[4], 0)

  // the smallest examples, that make sense for cross-linking
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, AASequence::fromString("HA"), 0, 2000.0, true, 1, 1);
  TEST_EQUAL(spec.size(), 1)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, AASequence::fromString("HA"), 1, 2000.0, true, 1, 1);
  TEST_EQUAL(spec.size(), 1)

  // loop link
  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, AASequence::fromString("PEPTIDESAREWEIRD"), 1, 2000.0, true, 1, 1, 14);
  TEST_EQUAL(spec.size(), 2)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, AASequence::fromString("PEPTIDESAREWEIRD"), 2, 2000.0, false, 1, 1, 14);
  TEST_EQUAL(spec.size(), 3)

  spec.clear(true);
  ptr->getXLinkIonSpectrum(spec, AASequence::fromString("PEPTIDESAREWEIRD"), 2, 2000.0, false, 1, 1, 13);
  TEST_EQUAL(spec.size(), 4)

  // test isotopic peaks
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 1); // not very useful combination, but it should at least run
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "false");
  param.setValue("add_metainfo", "false");
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);
  // 6 ion types with 4 charges each are expected
  TEST_EQUAL(spec.size(), 24)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 2); //
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);
  // 6 ion types with 4 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 48)

  spec.clear(true);
  param.setValue("add_isotopes", "true");
  param.setValue("max_isotope", 3); // not supported yet, but it should at least run (with the maximal possible number of peaks)
  ptr->setParameters(param);
  ptr->getXLinkIonSpectrum(spec, peptide, 3, 2000.0, true, 2, 5);
  // 6 ion types with 4 charges each are expected, each with a second isotopic peak
  TEST_EQUAL(spec.size(), 48)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  // for quick benchmarking of implementation chances
//  param = ptr->getParameters();
//  param.setValue("add_a_ions", "true");
//  param.setValue("add_b_ions", "true");
//  param.setValue("add_c_ions", "true");
//  param.setValue("add_x_ions", "true");
//  param.setValue("add_y_ions", "true");
//  param.setValue("add_z_ions", "true");
//  param.setValue("add_metainfo", "true");
//  ptr->setParameters(param);
//  AASequence tmp_peptide = AASequence::fromString("PEPTIDEPEPTIDEPEPTIDE");
//  for (Size i = 0; i != 1e4; ++i)
//  {
//    PeakSpectrum spec;
//    ptr->getXLinkIonSpectrum(spec, tmp_peptide, 9, 2000.0, false, 2, 5);
//  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

