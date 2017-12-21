
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
// $Maintainer: Timo Sachsenberg, Eugen Netz $
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
#include <OpenMS/CONCEPT/Constants.h>

///////////////////////////

START_TEST(TheoreticalSpectrumGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TheoreticalSpectrumGenerator* ptr = nullptr;
TheoreticalSpectrumGenerator* nullPointer = nullptr;

START_SECTION(TheoreticalSpectrumGenerator())
  ptr = new TheoreticalSpectrumGenerator();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source))
  TheoreticalSpectrumGenerator copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~TheoreticalSpectrumGenerator())
  delete ptr;
END_SECTION

ptr = new TheoreticalSpectrumGenerator();
AASequence peptide = AASequence::fromString("IFSQVGK");

START_SECTION(TheoreticalSpectrumGenerator& operator = (const TheoreticalSpectrumGenerator& tsg))
  TheoreticalSpectrumGenerator copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void getSpectrum(PeakSpectrum& spec, const AASequence& peptide, Int min_charge = 1, Int max_charge = 1))
  PeakSpectrum spec;
  ptr->getSpectrum(spec, peptide, 1, 1);
  TEST_EQUAL(spec.size(), 11)

  TOLERANCE_ABSOLUTE(0.001)

  double result[] = {/*114.091,*/ 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  ptr->getSpectrum(spec, peptide, 1, 2);
  TEST_EQUAL(spec.size(), 22)

  spec.clear(true);
  Param param(ptr->getParameters());
  param.setValue("add_first_prefix_ion", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, peptide, 1, 1);
  TEST_EQUAL(spec.size(), 12)

  double result2[] = {114.091, 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result2[i])
  }


  AASequence new_peptide = AASequence::fromString("DFPLANGER");
  /**  From http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html
   Seq    #       A            B            C            X            Y            Z         # (+1)
    D     1     88.03990    116.03481    133.06136       -        1018.49583   1001.46928    9
    F     2    235.10831    263.10323    280.12978    929.44815    903.46888    886.44233    8
    P     3    332.16108    360.15599    377.18254    782.37973    756.40047    739.37392    7
    I     4    445.24514    473.24005    490.26660    685.32697    659.34771    642.32116    6
    A     5    516.28225    544.27717    561.30372    572.24291    546.26364    529.23709    5
    N     6    630.32518    658.32009    675.34664    501.20579    475.22653    458.19998    4
    G     7    687.34664    715.34156    732.36811    387.16287    361.18360    344.15705    3
    E     8    816.38924    844.38415    861.41070    330.14140    304.16214    287.13559    2
    R     9    972.49035   1000.48526       -         201.09881    175.11955    158.09300    1
  **/
  double result_all[52-1] = {
   88.03990, 235.10831, 332.16108, 445.24514, 516.28225, 630.32518, 687.34664, 816.38924, /*972.49035, because TSG does not do A-ions of the full peptide*/
   116.03481,  263.10323,  360.15599,  473.24005,  544.27717,  658.32009,  715.34156,  844.38415, 1000.48526,
   133.06136, 280.12978, 377.18254, 490.26660, 561.30372, 675.34664, 732.36811, 861.41070,
   929.44815, 782.37973, 685.32697, 572.24291, 501.20579, 387.16287, 330.14140, 201.09881,
   1018.49583,  903.46888,  756.40047,  659.34771,  546.26364,  475.22653,  361.18360,  304.16214,  175.11955,
   1001.46928,  886.44233,  739.37392,  642.32116,  529.23709,  458.19998,  344.15705,  287.13559,  158.09300
  };
  std::sort(result_all,result_all+52-1);
  spec.clear(true);

  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_precursor_peaks", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, new_peptide, 1, 1);
  TEST_EQUAL(spec.size(), 52-1)

  vector<double> generated;
  for (Size i = 0; i != spec.size(); ++i)
  {
    generated.push_back(spec[i].getPosition()[0]);
  }

  std::sort(generated.begin(),generated.end());
  for (Size i = 0; i != generated.size(); ++i)
  {
    TEST_REAL_SIMILAR(generated[i], result_all[i])
  }

  // test loss creation and annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_precursor_peaks", "true");
  param.setValue("add_metainfo", "true");
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, peptide, 1, 1);
  TEST_EQUAL(spec.size(), 30)
  set<String> ion_names;
  // ions without losses
  ion_names.insert("b1+");
  ion_names.insert("x1+");
  ion_names.insert("b2+");
  ion_names.insert("x2+");
  ion_names.insert("b3+");
  ion_names.insert("x3+");
  ion_names.insert("b4+");
  ion_names.insert("x4+");
  ion_names.insert("b5+");
  ion_names.insert("x5+");
  ion_names.insert("b6+");
  ion_names.insert("x6+");

  // currently losses are generated independent of ion ladder type (b,y,...)
  // if an amino acid with potential loss is present in the prefix/suffix, then the loss is applied
  // if multiple amino acids with the same e.g. water loss are present in the prefix/suffix ion then the loss is only applied once
  ion_names.insert("x1-H3N1+");
  ion_names.insert("x2-H3N1+");
  ion_names.insert("x3-H3N1+");
  ion_names.insert("b3-H2O1+");
  ion_names.insert("x4-H3N1+");
  ion_names.insert("b4-H2O1+");
  ion_names.insert("b4-H3N1+");
  ion_names.insert("x5-H2O1+");
  ion_names.insert("x5-H3N1+");
  ion_names.insert("b5-H2O1+");
  ion_names.insert("b5-H3N1+");
  ion_names.insert("b6-H2O1+");
  ion_names.insert("b6-H3N1+");
  ion_names.insert("x6-H2O1+");
  ion_names.insert("x6-H3N1+");

  // precursors
  ion_names.insert("[M+H]-H2O+");
  ion_names.insert("[M+H]-NH3+");
  ion_names.insert("[M+H]+");

  PeakSpectrum::StringDataArray string_array = spec.getStringDataArrays().at(0);

  // check if all losses have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    String name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  // test for charges stored in IntegerDataArray
  PeakSpectrum charge3_spec;
  ptr->getSpectrum(charge3_spec, peptide, 1, 3);
  PeakSpectrum::IntegerDataArray charge_array = charge3_spec.getIntegerDataArrays().at(0);

  int charge_counts[3] = {0, 0, 0};
  for (Size i = 0; i != charge3_spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 27)
  TEST_EQUAL(charge_counts[1], 27)
  TEST_EQUAL(charge_counts[2], 30) // 3 more for [M+H], [M+H]-H20, [M+H]-NH3

  // test getSpectrum with one specific charge != 1
  spec.clear(true);
  ptr->getSpectrum(spec, peptide, 3, 3);
  TEST_EQUAL(spec.size(), 30)

  ion_names.clear();
  // ions without losses
  ion_names.insert("b1+++");
  ion_names.insert("x1+++");
  ion_names.insert("b2+++");
  ion_names.insert("x2+++");
  ion_names.insert("b3+++");
  ion_names.insert("x3+++");
  ion_names.insert("b4+++");
  ion_names.insert("x4+++");
  ion_names.insert("b5+++");
  ion_names.insert("x5+++");
  ion_names.insert("b6+++");
  ion_names.insert("x6+++");

  // currently losses are generated independent of ion ladder type (b,y,...)
  // if an amino acid with potential loss is present in the prefix/suffix, then the loss is applied
  // if multiple amino acids with the same e.g. water loss are present in the prefix/suffix ion then the loss is only applied once
  ion_names.insert("x1-H3N1+++");
  ion_names.insert("x2-H3N1+++");
  ion_names.insert("x3-H3N1+++");
  ion_names.insert("b3-H2O1+++");
  ion_names.insert("x4-H3N1+++");
  ion_names.insert("b4-H2O1+++");
  ion_names.insert("b4-H3N1+++");
  ion_names.insert("x5-H2O1+++");
  ion_names.insert("x5-H3N1+++");
  ion_names.insert("b5-H2O1+++");
  ion_names.insert("b5-H3N1+++");
  ion_names.insert("b6-H2O1+++");
  ion_names.insert("b6-H3N1+++");
  ion_names.insert("x6-H2O1+++");
  ion_names.insert("x6-H3N1+++");

  // precursors
  ion_names.insert("[M+H]-H2O+++");
  ion_names.insert("[M+H]-NH3+++");
  ion_names.insert("[M+H]+++");

  string_array = spec.getStringDataArrays().at(0);

  // check if all losses have been annotated
  for (Size i = 0; i != spec.size(); ++i)
  {
    String name = string_array[i];
    TEST_EQUAL(ion_names.find(name) != ion_names.end(), true)
  }

  charge_array = spec.getIntegerDataArrays().at(0);

  charge_counts[0] = 0;
  charge_counts[1] = 0;
  charge_counts[2] = 0;
  for (Size i = 0; i != spec.size(); ++i)
  {
    charge_counts[charge_array[i]-1]++;
  }
  TEST_EQUAL(charge_counts[0], 0)
  TEST_EQUAL(charge_counts[1], 0)
  TEST_EQUAL(charge_counts[2], 30)

  // AbundantImmoniumIons test
  param = ptr->getParameters();
  param.setValue("add_b_ions", "false");
  param.setValue("add_x_ions", "false");
  param.setValue("add_precursor_peaks", "false");
  param.setValue("add_metainfo", "false");
  param.setValue("add_losses", "false");
  param.setValue("add_abundant_immonium_ions", "true");
  ptr->setParameters(param);
  spec.clear(true);
  ptr->getSpectrum(spec, AASequence::fromString("HFYLWCP"), 1, 1);
  TEST_EQUAL(spec.size(), 7)
  TEST_REAL_SIMILAR(spec[0].getPosition()[0], 70.0656)
  TEST_REAL_SIMILAR(spec[1].getPosition()[0], 76.0221)
  TEST_REAL_SIMILAR(spec[2].getPosition()[0], 86.09698)
  TEST_REAL_SIMILAR(spec[3].getPosition()[0], 110.0718)
  TEST_REAL_SIMILAR(spec[4].getPosition()[0], 120.0813)
  TEST_REAL_SIMILAR(spec[5].getPosition()[0], 136.0762)
  TEST_REAL_SIMILAR(spec[6].getPosition()[0], 159.0922)

  spec.clear(true);
  ptr->getSpectrum(spec, AASequence::fromString("H"), 1, 1);
  TEST_EQUAL(spec.size(), 1)

  spec.clear(true);
  ptr->getSpectrum(spec, AASequence::fromString("A"), 1, 1);
  TEST_EQUAL(spec.size(), 0)


//  // for quick benchmarking of implementation chances
//  param = ptr->getParameters();
//  param.setValue("add_first_prefix_ion", "true");
//  param.setValue("add_a_ions", "true");
//  param.setValue("add_b_ions", "true");
//  param.setValue("add_c_ions", "true");
//  param.setValue("add_x_ions", "true");
//  param.setValue("add_y_ions", "true");
//  param.setValue("add_z_ions", "true");
//  param.setValue("add_precursor_peaks", "true");
//  param.setValue("add_metainfo", "true");
//  param.setValue("add_losses", "true");
//  ptr->setParameters(param);
//  peptide = AASequence::fromString("PEPTIDEPEPTIDEPEPTIDE");
//  for (Size i = 0; i != 1e4; ++i)
//  {
//    PeakSpectrum spec;
//    ptr->getSpectrum(spec, peptide, 1, 3);
//  }

END_SECTION

START_SECTION(([EXTRA] bugfix test where losses lead to formulae with negative element frequencies))
{
  AASequence tmp_aa = AASequence::fromString("RDAGGPALKK");
  PeakSpectrum tmp;
  TheoreticalSpectrumGenerator t_gen;
  Param params;

  params.setValue("add_isotopes", "true");
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_a_ions", "true");
  t_gen.setParameters(params);

  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp.size(), 212)
}
END_SECTION

START_SECTION(([EXTRA] test monomer extreme case))
{
  AASequence tmp_aa = AASequence::fromString("R");
  PeakSpectrum tmp;
  TheoreticalSpectrumGenerator t_gen;
  Param params;

  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_x_ions", "true");
  t_gen.setParameters(params);
  TEST_EXCEPTION(Exception::InvalidSize, t_gen.getSpectrum(tmp, tmp_aa, 1, 1));

  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_x_ions", "false");
  params.setValue("add_c_ions", "true");
  t_gen.setParameters(params);
  TEST_EXCEPTION(Exception::InvalidSize, t_gen.getSpectrum(tmp, tmp_aa, 1, 1));

  params.setValue("add_x_ions", "false");
  params.setValue("add_c_ions", "false");
  params.setValue("add_precursor_peaks", "true");
  t_gen.setParameters(params);
  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp.size(), 3)
}
END_SECTION

START_SECTION(([EXTRA] test isotope clusters for all peak types))
{
  AASequence tmp_aa = AASequence::fromString("ARRGH");
  PeakSpectrum spec;
  TheoreticalSpectrumGenerator t_gen;
  Param params;
  params.setValue("add_isotopes", "true");
  params.setValue("max_isotope", 2);
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);

  // isotope cluster for y-ions
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 8)

  TOLERANCE_ABSOLUTE(0.001)
  double neutron_shift = Constants::NEUTRON_MASS_U;

  // 4 monoisotopic masses, 4 second peaks with added neutron mass / 2
  double result[] = {78.54206, 107.05279, 185.10335, 263.15390, 78.54206+(neutron_shift/2), 107.05279+(neutron_shift/2), 185.10335+(neutron_shift/2), 263.15390+(neutron_shift/2)};
  std::sort(result, result+8);
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  // isotope cluster for losses
  spec.clear(true);
  params.setValue("add_losses", "true");
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);
  t_gen.getSpectrum(spec, tmp_aa, 1, 2);
  TEST_EQUAL(spec.size(), 40)

  double proton_shift = Constants::PROTON_MASS_U;
  // 10 monoisotopic peaks with charge=1, 10 second peaks, 20 with charge=2
  double result_losses[] = {156.07675, 213.09821, 325.18569, 327.17753, 352.17278, 369.19932, 481.28680, 483.27864, 508.27389, 525.30044,
                                           156.07675+neutron_shift, 213.09821+neutron_shift, 325.18569+neutron_shift, 327.17753+neutron_shift, 352.17278+neutron_shift, 369.19932+neutron_shift, 481.28680+neutron_shift, 483.27864+neutron_shift, 508.27389+neutron_shift, 525.30044+neutron_shift,
                                           (156.07675+proton_shift)/2, (213.09821+proton_shift)/2, (325.18569+proton_shift)/2, (327.17753+proton_shift)/2, (352.17278+proton_shift)/2, (369.19932+proton_shift)/2, (481.28680+proton_shift)/2, (483.27864+proton_shift)/2, (508.27389+proton_shift)/2, (525.30044+proton_shift)/2,
                                           (156.07675+proton_shift)/2+(neutron_shift/2), (213.09821+proton_shift)/2+(neutron_shift/2), (325.18569+proton_shift)/2+(neutron_shift/2), (327.17753+proton_shift)/2+(neutron_shift/2), (352.17278+proton_shift)/2+(neutron_shift/2),
                                           (369.19932+proton_shift)/2+(neutron_shift/2), (481.28680+proton_shift)/2+(neutron_shift/2), (483.27864+proton_shift)/2+(neutron_shift/2), (508.27389+proton_shift)/2+(neutron_shift/2), (525.30044+proton_shift)/2+(neutron_shift/2)};
  std::sort(result_losses, result_losses+40);
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result_losses[i])
  }

  // isotope cluster for precurser peaks with losses
  spec.clear(true);
  params.setValue("add_precursor_peaks", "true");
  params.setValue("add_b_ions", "false");
  params.setValue("add_y_ions", "false");

  t_gen.setParameters(params);
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 6)

  // 3 monoisitopic peaks, 3 second peaks
  double result_precursors[] = {(578.32698+proton_shift)/2, (579.31100+proton_shift)/2, (596.33755+proton_shift)/2,
                                                  (578.32698+proton_shift)/2+(neutron_shift/2), (579.31100+proton_shift)/2+(neutron_shift/2), (596.33755+proton_shift)/2+(neutron_shift/2)};
  std::sort(result_precursors, result_precursors+6);
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result_precursors[i])
  }


}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
