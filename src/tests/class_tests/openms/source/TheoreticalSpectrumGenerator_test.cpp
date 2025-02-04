// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
  delete ptr;
END_SECTION

START_SECTION(TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source))
  ptr = new TheoreticalSpectrumGenerator();
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

  /**  From http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html
     Fragment Ion Table, monoisotopic masses

     Seq    #       A            B            C            X            Y            Z         # (+1)

     I     1     86.09647    114.09139    131.11793       -         778.44581    761.42036    7
     F     2    233.16488    261.15980    278.18635    691.34101    665.36174    648.33629    6
     S     3    320.19691    348.19183    365.21838    544.27260    518.29333    501.26788    5
     Q     4    448.25549    476.25040    493.27695    457.24057    431.26130    414.23585    4
     V     5    547.32390    575.31882    592.34537    329.18199    303.20273    286.17727    3
     G     6    604.34537    632.34028    649.36683    230.11358    204.13431    187.10886    2
     K     7    732.44033    760.43524       -         173.09211    147.11285    130.08740    1

  **/
  double result[] = {/*114.091,*/ 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
  std::vector<double> result_x = { 691.34101, 544.27260, 457.24057, 329.18199, 230.11358, 173.09211 };
  std::vector<double> result_x_losses = {
      691.34101 - 17.026549095700005,
      691.34101 - 18.01056506379996,
      691.34101,
      544.27260 - 17.026549095700005,
      544.27260 - 18.01056506379996,
      544.27260,
      457.24057 - 17.026549095700005,
      457.24057,
      329.18199 - 17.026549095700005,
      329.18199,
      230.11358 - 17.026549095700005,
      230.11358,
      173.09211 - 17.026549095700005,
      173.09211 };
  std::sort(result_x.begin(), result_x.end());
  std::sort(result_x_losses.begin(), result_x_losses.end());

  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  TEST_EQUAL(spec.getMSLevel(), 2);
  TEST_EQUAL(spec.getType(), MSSpectrum::SpectrumSettings::CENTROID);
  TEST_REAL_SIMILAR(peptide.getMZ(2, Residue::Full), spec.getPrecursors()[0].getMZ());

  spec.clear(true);
  ptr->getSpectrum(spec, peptide, 1, 2);
  TEST_EQUAL(spec.size(), 22)

  TEST_REAL_SIMILAR(peptide.getMZ(3, Residue::Full), spec.getPrecursors()[0].getMZ());

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

  std::vector<double> result_bx = {
   116.03481,  263.10323,  360.15599,  473.24005,  544.27717,  658.32009,  715.34156,  844.38415, 1000.48526,
   929.44815, 782.37973, 685.32697, 572.24291, 501.20579, 387.16287, 330.14140, 201.09881,
  };
  std::sort(result_bx.begin(),result_bx.end());

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

  TEST_REAL_SIMILAR(new_peptide.getMZ(2, Residue::Full), spec.getPrecursors()[0].getMZ());

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
  param.setValue("add_b_ions", "false");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_precursor_peaks", "false");
  param.setValue("add_metainfo", "false");
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, peptide, 1, 1);
  TEST_EQUAL(spec.size(), 14)

  generated.clear();
  for (Size i = 0; i != spec.size(); ++i)
  {
    generated.push_back(spec[i].getPosition()[0]);
  }
  for (Size i = 0; i != generated.size(); ++i)
  {
    TEST_REAL_SIMILAR(generated[i], result_x_losses[i])
  }

  // test loss creation and annotation
  spec.clear(true);
  param = ptr->getParameters();
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "false");
  param.setValue("add_c_ions", "false");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "false");
  param.setValue("add_precursor_peaks", "false");
  param.setValue("add_metainfo", "true");
  param.setValue("add_losses", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spec, peptide, 1, 1);
  TEST_EQUAL(spec.size(), 14)

  generated.clear();
  for (Size i = 0; i != spec.size(); ++i)
  {
    generated.push_back(spec[i].getPosition()[0]);
  }
  for (Size i = 0; i != generated.size(); ++i)
  {
    TEST_REAL_SIMILAR(generated[i], result_x_losses[i])
  }

  std::sort(generated.begin(),generated.end());

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

  TEST_REAL_SIMILAR(peptide.getMZ(4, Residue::Full), spec.getPrecursors()[0].getMZ());

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

  spec.clear(true);
  ptr->getSpectrum(spec, peptide, 1, 1, 4);
  ptr->getSpectrum(spec, new_peptide, 1, 3);
  ABORT_IF(spec.getPrecursors().size() != 2);
  TEST_REAL_SIMILAR(spec.getPrecursors()[0].getMZ(), peptide.getMZ(4));
  TEST_EQUAL(spec.getPrecursors()[0].getCharge(), 4);
  TEST_REAL_SIMILAR(spec.getPrecursors()[1].getMZ(), new_peptide.getMZ(4));
  TEST_EQUAL(spec.getPrecursors()[1].getCharge(), 4);

  spec.clear(true);

  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ptr->getSpectrum(spec, peptide, 1, 2, 1), "'precursor_charge' has to be higher than or equal to 'max_charge'.");

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

START_SECTION(static MSSpectrum generateSpectrum(const Precursor::ActivationMethod& fm, const AASequence& seq, int precursor_charge))
  MSSpectrum spec;
  Precursor prec;

  // Test CID/HCID
  spec = TheoreticalSpectrumGenerator::generateSpectrum(Precursor::ActivationMethod::CID, AASequence::fromString("HFYLWCP"), 1);
  ABORT_IF(spec.size() != 11);
  TEST_REAL_SIMILAR(spec[0].getPosition()[0], 116.0706);
  TEST_REAL_SIMILAR(spec[1].getPosition()[0], 219.0797);
  TEST_REAL_SIMILAR(spec[2].getPosition()[0], 285.1346);
  TEST_REAL_SIMILAR(spec[3].getPosition()[0], 405.1591);
  TEST_REAL_SIMILAR(spec[4].getPosition()[0], 448.1979);
  TEST_REAL_SIMILAR(spec[5].getPosition()[0], 518.2431);
  TEST_REAL_SIMILAR(spec[6].getPosition()[0], 561.2819);
  TEST_REAL_SIMILAR(spec[7].getPosition()[0], 681.3064);
  TEST_REAL_SIMILAR(spec[8].getPosition()[0], 747.3613);
  TEST_REAL_SIMILAR(spec[9].getPosition()[0], 828.3749);
  TEST_REAL_SIMILAR(spec[10].getPosition()[0], 850.3704);

  spec.clear(true);

  // Test ECD/ETD
  spec = TheoreticalSpectrumGenerator::generateSpectrum(Precursor::ActivationMethod::ECD, AASequence::fromString("HFYLWCP"), 1);
  TEST_EQUAL(spec.size(), 17);

  TEST_REAL_SIMILAR(spec[0].getPosition()[0], 100.0518816);
  TEST_REAL_SIMILAR(spec[1].getPosition()[0], 101.0597067);
  TEST_REAL_SIMILAR(spec[2].getPosition()[0], 203.0610665);
  TEST_REAL_SIMILAR(spec[3].getPosition()[0], 204.0688916);
  TEST_REAL_SIMILAR(spec[4].getPosition()[0], 302.1611520);
  TEST_REAL_SIMILAR(spec[5].getPosition()[0], 389.1403798);
  TEST_REAL_SIMILAR(spec[6].getPosition()[0], 390.1482049);
  TEST_REAL_SIMILAR(spec[7].getPosition()[0], 465.2244813);
  TEST_REAL_SIMILAR(spec[8].getPosition()[0], 502.2244442);
  TEST_REAL_SIMILAR(spec[9].getPosition()[0], 503.2322692);
  TEST_REAL_SIMILAR(spec[10].getPosition()[0], 578.3085457);
  // ...

  spec.clear(true);

  // Test precursor_charge > 2
  spec = TheoreticalSpectrumGenerator::generateSpectrum(Precursor::ActivationMethod::HCID, AASequence::fromString("PEP"), 3);
  TEST_EQUAL(spec.size(), 8);
  TEST_REAL_SIMILAR(spec[0].getPosition()[0], 58.5389);
  TEST_REAL_SIMILAR(spec[1].getPosition()[0], 100.0574);
  TEST_REAL_SIMILAR(spec[2].getPosition()[0], 114.0549);
  TEST_REAL_SIMILAR(spec[3].getPosition()[0], 116.0706);
  TEST_REAL_SIMILAR(spec[4].getPosition()[0], 123.0602);
  TEST_REAL_SIMILAR(spec[5].getPosition()[0], 199.1077);
  TEST_REAL_SIMILAR(spec[6].getPosition()[0], 227.1026);
  TEST_REAL_SIMILAR(spec[7].getPosition()[0], 245.1131);

  // Test not supported activation method
  TEST_EXCEPTION(Exception::InvalidParameter, TheoreticalSpectrumGenerator::generateSpectrum(Precursor::ActivationMethod::SORI, AASequence::fromString("PEP"), 1));

END_SECTION

START_SECTION(([EXTRA] bugfix test where losses lead to formulae with negative element frequencies))
{
  // this tests for the loss of CONH2 on Arginine, however it is not clear how
  // this loss would occur in the first place.
  AASequence tmp_aa = AASequence::fromString("RDAGGPALKK");
  PeakSpectrum tmp;
  TheoreticalSpectrumGenerator t_gen;
  Param params;

  params.setValue("isotope_model", "coarse");
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_a_ions", "true");
  t_gen.setParameters(params);

  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp.size(), 212)

  tmp.clear(true);
  params.setValue("isotope_model", "coarse");
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "false");
  params.setValue("add_a_ions", "true");
  t_gen.setParameters(params);
  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp_aa[0].hasNeutralLoss(), true)
  TEST_EQUAL(tmp.size(), 198)

  tmp_aa = AASequence::fromString("RDK");
  tmp.clear(true);
  params.setValue("isotope_model", "none");
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_a_ions", "true");
  params.setValue("add_b_ions", "false");
  params.setValue("add_y_ions", "false");
  params.setValue("add_metainfo", "true");
  t_gen.setParameters(params);

  TEST_EQUAL(tmp.size(), 0)
  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  // TEST_EQUAL(tmp.size(), 8)

  tmp.clear(true);
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_a_ions", "true");
  params.setValue("add_b_ions", "false");
  params.setValue("add_y_ions", "false");
  params.setValue("add_metainfo", "false");
  t_gen.setParameters(params);

  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  // TEST_EQUAL(tmp.size(), 8)
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
  params.setValue("isotope_model", "coarse");
  params.setValue("max_isotope", 2);
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);

  // isotope cluster for y-ions
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 8)

  TOLERANCE_ABSOLUTE(0.001)
  double neutron_shift = Constants::C13C12_MASSDIFF_U;

  // 4 monoisotopic masses, 4 second peaks with added neutron mass / 2
  std::vector<double> result = {
    78.54206,
    107.05279,
    185.10335,
    263.15390,

    78.54206+(neutron_shift/2),
    107.05279+(neutron_shift/2),
    185.10335+(neutron_shift/2),
    263.15390+(neutron_shift/2)
  };

  std::sort(result.begin(), result.end());
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  params.setValue("isotope_model", "fine");
  params.setValue("max_isotope", 2);
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);

  // isotope cluster for y-ions
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 10)

  result = {
    78.54206,
    107.05279,
    185.10335,
    263.15390,

    // 405: POS: 78.54256367545 INT: 0.921514272689819
    79.04424117545,
    // 405: POS: 107.0532957233 INT: 0.89608770608902
    107.5549732233,
    // 405: POS: 185.1038514147 INT: 0.824628114700317
    185.6023689147,
    185.6055289147,
    // 405: POS: 263.1544071061 INT: 0.758867204189301
    263.6529246061,
    263.6560846061,
  };
  std::sort(result.begin(), result.end());
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  params.setValue("isotope_model", "fine");
  params.setValue("max_isotope", 2);
  params.setValue("max_isotope_probability", 0.20);
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);

  // isotope cluster for y-ions
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 5)

  result = {
    78.54206,
    107.05279,
    185.10335,
    263.15390,

    // 405: POS: 78.54256367545 INT: 0.921514272689819
    // 405: POS: 107.0532957233 INT: 0.89608770608902
    // 405: POS: 185.1038514147 INT: 0.824628114700317
    // 405: POS: 263.1544071061 INT: 0.758867204189301
    263.6560846061,
  };
  std::sort(result.begin(), result.end());
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
  }

  spec.clear(true);
  params.setValue("isotope_model", "fine");
  params.setValue("max_isotope", 2);
  params.setValue("max_isotope_probability", 0.01);
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);

  // isotope cluster for y-ions
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  // TEST_EQUAL(spec.size(), 34)

  // isotope cluster for losses
  spec.clear(true);
  params.setValue("isotope_model", "coarse");
  params.setValue("add_losses", "true");
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);
  t_gen.getSpectrum(spec, tmp_aa, 1, 2);
  TEST_EQUAL(spec.size(), 40)

  double proton_shift = Constants::PROTON_MASS_U;
  // 10 monoisotopic peaks with charge=1, 10 second peaks, 20 with charge=2
  std::vector<double> result_losses = { 156.07675, 213.09821, 325.18569, 327.17753, 352.17278, 369.19932, 481.28680, 483.27864, 508.27389, 525.30044,
	   156.07675+neutron_shift, 213.09821+neutron_shift, 325.18569+neutron_shift, 327.17753+neutron_shift, 352.17278+neutron_shift, 369.19932+neutron_shift, 481.28680+neutron_shift, 483.27864+neutron_shift, 508.27389+neutron_shift, 525.30044+neutron_shift,
	  (156.07675+proton_shift)/2, (213.09821+proton_shift)/2, (325.18569+proton_shift)/2, (327.17753+proton_shift)/2, (352.17278+proton_shift)/2, (369.19932+proton_shift)/2, (481.28680+proton_shift)/2, (483.27864+proton_shift)/2, (508.27389+proton_shift)/2, (525.30044+proton_shift)/2,
	  (156.07675+proton_shift)/2+(neutron_shift/2), (213.09821+proton_shift)/2+(neutron_shift/2), (325.18569+proton_shift)/2+(neutron_shift/2), (327.17753+proton_shift)/2+(neutron_shift/2), (352.17278+proton_shift)/2+(neutron_shift/2),
	  (369.19932+proton_shift)/2+(neutron_shift/2), (481.28680+proton_shift)/2+(neutron_shift/2), (483.27864+proton_shift)/2+(neutron_shift/2), (508.27389+proton_shift)/2+(neutron_shift/2), (525.30044+proton_shift)/2+(neutron_shift/2)};

  std::sort(result_losses.begin(), result_losses.end());
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result_losses[i])
  }
  result_losses = { 0.927642, 0.0723581}; // check intensity
  for (Size i = 0; i != 2; ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getIntensity(), result_losses[i])
  }

  // last two entries:
  TEST_REAL_SIMILAR( spec[ spec.size() -2 ].getMZ(), 525.30044)
  TEST_REAL_SIMILAR( spec[ spec.size() -1 ].getMZ(), 526.304)

  spec.clear(true);
  params.setValue("isotope_model", "fine");
  params.setValue("max_isotope_probability", 0.05);
  params.setValue("add_losses", "true");
  params.setValue("add_b_ions", "false");
  t_gen.setParameters(params);
  t_gen.getSpectrum(spec, tmp_aa, 1, 2);
  TEST_EQUAL(spec.size(), 50)

  result_losses = { 78.5426, 79.0442, 107.0532, 107.5549};
  for (Size i = 0; i != 4; ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getMZ(), result_losses[i])
  }
  result_losses = { 0.921514, 0.0598011, 0.896088, 0.0775347}; // check intensity
  for (Size i = 0; i != 4; ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getIntensity(), result_losses[i])
  }

  // last entries
  TEST_REAL_SIMILAR( spec[ spec.size() -5 ].getMZ(), 509.271)
  TEST_REAL_SIMILAR( spec[ spec.size() -4 ].getMZ(), 509.277)
  TEST_REAL_SIMILAR( spec[ spec.size() -3 ].getMZ(), 525.301)
  TEST_REAL_SIMILAR( spec[ spec.size() -2 ].getMZ(), 526.298)
  TEST_REAL_SIMILAR( spec[ spec.size() -1 ].getMZ(), 526.304)

  // isotope cluster for precursor peaks with losses
  spec.clear(true);
  params.setValue("add_precursor_peaks", "true");
  params.setValue("isotope_model", "coarse");
  params.setValue("add_b_ions", "false");
  params.setValue("add_y_ions", "false");

  t_gen.setParameters(params);
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 6)

  // 3 monoisotopic peaks, 3 second peaks
  double result_precursors[] = {
    (578.32698+proton_shift)/2,
      (579.31100+proton_shift)/2,
        (596.33755+proton_shift)/2,

    (578.32698+proton_shift)/2+(neutron_shift/2),
      (579.31100+proton_shift)/2+(neutron_shift/2),
        (596.33755+proton_shift)/2+(neutron_shift/2)};

  std::sort(result_precursors, result_precursors+6);
  for (Size i = 0; i != spec.size(); ++i)
  {
    TEST_REAL_SIMILAR(spec[i].getPosition()[0], result_precursors[i])
  }

  spec.clear(true);
  params.setValue("add_precursor_peaks", "true");
  params.setValue("isotope_model", "fine");
  params.setValue("add_b_ions", "false");
  params.setValue("add_y_ions", "false");

  t_gen.setParameters(params);
  t_gen.getSpectrum(spec, tmp_aa, 2, 2);
  TEST_EQUAL(spec.size(), 12)

  TEST_REAL_SIMILAR(spec[0].getMZ(), (578.32698+proton_shift)/2 )
  TEST_REAL_SIMILAR(spec[1].getMZ(), (579.31100+proton_shift)/2 )
  TEST_REAL_SIMILAR(spec[11].getMZ(), (598.34481333943+proton_shift)/2 )
}
END_SECTION

START_SECTION(([EXTRA] test SpectrumAnnotator ))
{
  // use same params as SpectrumAnnotator
  AASequence tmp_aa = AASequence::fromString("IALSRPNVEVVALNDPFITNDYAAYM(Oxidation)FK");
  PeakSpectrum tmp;
  TheoreticalSpectrumGenerator t_gen;
  Param tgp;
  tgp.setValue("add_metainfo", "true");
  tgp.setValue("add_losses", "true");
  tgp.setValue("add_precursor_peaks", "true");
  tgp.setValue("add_abundant_immonium_ions", "true");
  tgp.setValue("add_first_prefix_ion", "true");
  tgp.setValue("add_y_ions", "true");
  tgp.setValue("add_b_ions", "true");
  tgp.setValue("add_a_ions", "true");
  tgp.setValue("add_x_ions", "true");
  t_gen.setParameters(tgp);
  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp.size(), 465)

  tmp.clear(true);
  tgp.setValue("add_metainfo", "true");
  tgp.setValue("add_losses", "true");
  tgp.setValue("add_precursor_peaks", "false");
  tgp.setValue("add_abundant_immonium_ions", "false");
  tgp.setValue("add_first_prefix_ion", "false");
  tgp.setValue("add_y_ions", "false");
  tgp.setValue("add_b_ions", "false");
  tgp.setValue("add_a_ions", "true");
  tgp.setValue("add_x_ions", "false");
  t_gen.setParameters(tgp);
  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp.size(), 121)

  // for (Size k = 0; k < tmp.size(); k++)
  // {
  //   std::cout << tmp[k] << " -- " << tmp.getStringDataArrays()[0][k] << std::endl;
  // }

}
END_SECTION

START_SECTION(([EXTRA] test first prefix loss))
{
  AASequence tmp_aa = AASequence::fromString("RDAGGPALKK");
  PeakSpectrum tmp;
  TheoreticalSpectrumGenerator t_gen;
  Param params;

  params.setValue("isotope_model", "none");
  params.setValue("add_losses", "true");
  params.setValue("add_first_prefix_ion", "true");
  params.setValue("add_a_ions", "true");
  params.setValue("add_metainfo", "true");
  t_gen.setParameters(params);

  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp.size(), 107)

  auto anno = tmp.getStringDataArrays()[0];
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1+") != anno.end(), true)
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1-H3N1+") != anno.end(), true)
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1-C1H2N2+") != anno.end(), true)
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1-C1H2N1O1+") != anno.end(), true)

  // test without prefix ion (but still requires correct losses elsewhere)
  tmp.clear(true);
  params.setValue("add_first_prefix_ion", "false");
  t_gen.setParameters(params);
  t_gen.getSpectrum(tmp, tmp_aa, 1, 1);
  TEST_EQUAL(tmp_aa[0].hasNeutralLoss(), true)
  TEST_EQUAL(tmp.size(), 99) // missing a1 and b1 ions as well as their losses -H3N1+ C1H2N2+ -C1H2N1O1+

  anno = tmp.getStringDataArrays()[0];
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1+") == anno.end(), true)
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1-H3N1+") == anno.end(), true)
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1-C1H2N2+") == anno.end(), true)
  TEST_EQUAL(std::find(anno.begin(), anno.end(), "b1-C1H2N1O1+") == anno.end(), true)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
