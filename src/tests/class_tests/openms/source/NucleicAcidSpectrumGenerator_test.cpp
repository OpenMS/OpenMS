// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>

///////////////////////////

START_TEST(NucleicAcidSpectrumGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

NucleicAcidSpectrumGenerator* ptr = nullptr;
NucleicAcidSpectrumGenerator* null_ptr = nullptr;

START_SECTION(NucleicAcidSpectrumGenerator())
  ptr = new NucleicAcidSpectrumGenerator();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION(NucleicAcidSpectrumGenerator(const NucleicAcidSpectrumGenerator& source))
  NucleicAcidSpectrumGenerator copy(*ptr);
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(NucleicAcidSpectrumGenerator& operator=(const TheoreticalSpectrumGenerator& source))
  NucleicAcidSpectrumGenerator copy;
  copy = *ptr;
  TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~NucleicAcidSpectrumGenerator())
  delete ptr;
END_SECTION

ptr = new NucleicAcidSpectrumGenerator();

START_SECTION((void getSpectrum(MSSpectrum& spectrum, const NASequence& oligo, Int min_charge, Int max_charge) const))
{
  // fragment ion data from Ariadne (ariadne.riken.jp):
  NASequence seq = NASequence::fromString("[m1A]UCCACAGp");
  ABORT_IF(abs(seq.getMonoWeight() - 2585.3800) > 0.01);
  vector<double> aminusB_ions = {113.0244, 456.0926, 762.1179, 1067.1592,
                                 1372.2005, 1701.2530, 2006.2943, 2335.3468};
  vector<double> a_ions = {262.0946, 568.1199, 873.1612, 1178.2024, 1507.2550,
                           1812.2962, 2141.3488, 2486.3962};
  vector<double> b_ions = {280.1051, 586.1304, 891.1717, 1196.2130, 1525.2655,
                           1830.3068, 2159.3593, 2504.4068};
  vector<double> c_ions = {342.0609, 648.0862, 953.1275, 1258.1688, 1587.2213,
                           1892.2626, 2221.3151};
  vector<double> d_ions = {360.0715, 666.0968, 971.1380, 1276.1793, 1605.2319,
                           1910.2731, 2239.3257};
  vector<double> w_ions = {442.0171, 771.0696, 1076.1109, 1405.1634, 1710.2047,
                           2015.2460, 2321.2713};
  vector<double> x_ions = {424.0065, 753.0590, 1058.1003, 1387.1528, 1692.1941,
                           1997.2354, 2303.2607};
  vector<double> y_ions = {362.0507, 691.1032, 996.1445, 1325.1970, 1630.2383,
                           1935.2796, 2241.3049};
  vector<double> z_ions = {344.0402, 673.0927, 978.1340, 1307.1865, 1612.2278,
                           1917.2691, 2223.2944};

  Param param = ptr->getDefaults();
  param.setValue("add_metainfo", "true");
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_b_ions", "false");
  param.setValue("add_y_ions", "false");

  MSSpectrum spectrum;
  param.setValue("add_a-B_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), aminusB_ions.size() - 1); // last one is missing
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), aminusB_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_a-B_ions", "false");
  param.setValue("add_a_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), a_ions.size() - 1); // last one is missing
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), a_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), b_ions.size() - 1); // last one is missing
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), b_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_b_ions", "false");
  param.setValue("add_c_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), c_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), c_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_c_ions", "false");
  param.setValue("add_d_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), d_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), d_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_d_ions", "false");
  param.setValue("add_w_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), w_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), w_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_w_ions", "false");
  param.setValue("add_x_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), x_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), x_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), y_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), y_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), z_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), z_ions[i]);
  }

  seq = NASequence::fromString("[m1A]UCCACA[G*]p"); //Terminal thiol replacement shouldn't change any masses
  TEST_REAL_SIMILAR(seq.getMonoWeight(), 2585.3800);

  //repeat the above with internal Thiols

  seq = NASequence::fromString("[m1A]UC[C*]AC[A*]Gp"); //Terminal thiol replacement shouldn't change any masses
  ABORT_IF(abs(seq.getMonoWeight() - 2617.334342) > 0.01);
  aminusB_ions = {113.0244, 456.0926, 762.1179, 1067.1592,
                                 1388.1777, 1717.2302, 2022.27147, 2367.3011};
  a_ions = {262.0946, 568.1199, 873.1612, 1178.2024, 1523.2321,
                           1828.2733, 2157.3259, 2518.3505};
  b_ions = {280.1051, 586.1304, 891.1717, 1196.2130, 1541.2426,
                           1846.2839, 2175.3365, 2536.3611};
  c_ions = {342.0609, 648.0862, 953.1275, 1274.1458, 1603.1984,
                           1908.2397, 2253.2694};
  d_ions = {360.0715, 666.0968, 971.1380, 1292.1564, 1621.2090,
                           1926.2502, 2271.2800};
  w_ions = {457.9942, 787.0468, 1092.0881, 1437.1178, 1742.1591,
                           2047.2003, 2353.2256};
  x_ions = {439.9837, 769.0362, 1074.0775, 1419.1072, 1724.1485,
                           2029.1898, 2335.2150};
  y_ions = {362.0507, 707.0805, 1012.1217, 1341.1743, 1662.1927,
                           1967.2340, 2273.2593};
  z_ions = {344.0402, 689.0699, 994.1112, 1323.1637, 1644.1822,
                           1949.2234, 2255.2487};

  param = ptr->getDefaults();
  param.setValue("add_metainfo", "true");
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_b_ions", "false");
  param.setValue("add_y_ions", "false");

  spectrum.clear(true);

  param.setValue("add_a-B_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), aminusB_ions.size() - 1); // last one is missing
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), aminusB_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_a-B_ions", "false");
  param.setValue("add_a_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), a_ions.size() - 1); // last one is missing
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), a_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_a_ions", "false");
  param.setValue("add_b_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), b_ions.size() - 1); // last one is missing
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), b_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_b_ions", "false");
  param.setValue("add_c_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), c_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), c_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_c_ions", "false");
  param.setValue("add_d_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), d_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), d_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_d_ions", "false");
  param.setValue("add_w_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), w_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), w_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_w_ions", "false");
  param.setValue("add_x_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), x_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), x_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_x_ions", "false");
  param.setValue("add_y_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), y_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), y_ions[i]);
  }

  spectrum.clear(true);
  param.setValue("add_y_ions", "false");
  param.setValue("add_z_ions", "true");
  ptr->setParameters(param);
  ptr->getSpectrum(spectrum, seq, -1, -1);
  TEST_EQUAL(spectrum.size(), z_ions.size());
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getMZ(), z_ions[i]);
  }


}
END_SECTION


START_SECTION((void getMultipleSpectra(std::map<Int, MSSpectrum>& spectra, const NASequence& oligo, const std::set<Int>& charges, Int base_charge = 1) const))
{
  NucleicAcidSpectrumGenerator gen;
  Param param = gen.getParameters();
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_metainfo", "true");
  // param.setValue("add_precursor_peaks", "true"); // yes or no?
  param.setValue("add_a_ions", "true");
  param.setValue("add_b_ions", "true");
  param.setValue("add_c_ions", "true");
  param.setValue("add_d_ions", "true");
  param.setValue("add_w_ions", "true");
  param.setValue("add_x_ions", "true");
  param.setValue("add_y_ions", "true");
  param.setValue("add_z_ions", "true");
  param.setValue("add_a-B_ions", "true");

  NASequence seq = NASequence::fromString("[m1A]UCCACAGp");
  set<Int> charges = {-1, -3, -5};
  // get spectra individually:
  vector<MSSpectrum> compare(charges.size());
  Size index = 0;
  for (Int charge : charges)
  {
    gen.getSpectrum(compare[index], seq, -1, charge);
    index++;
  }
  // now all together:
  map<Int, MSSpectrum> spectra;
  gen.getMultipleSpectra(spectra, seq, charges, -1);
  // compare:
  TEST_EQUAL(compare.size(), spectra.size());
  index = 0;
  for (const auto& pair : spectra)
  {
    TEST_EQUAL(compare[index] == pair.second, true);
    index++;
  }
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
