// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/MISC/EmgGradientDescent.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

// TODO: remove helper function and its calls
void geogebra_print_execute(const double h, const double mu, const double sigma, const double tau)
{
  std::cout << "\nGEOGEBRA: Execute[{\"h = " << h << "\", \"mu = " << mu << "\",\"sigma = " << sigma << "\", \"tau = " << tau << "\"}]\n\n";
}

START_TEST(EmgGradientDescent, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EmgGradientDescent* ptr = nullptr;
EmgGradientDescent* null_ptr = nullptr;

// Toy chromatogram
// data is taken from raw LC-MS/MS data points acquired for L-Glutamate in RBCs
const vector<double> position = {
  2.23095,2.239716667,2.248866667,2.25765,2.266416667,
  2.275566667,2.2847,2.293833333,2.304066667,2.315033333,2.325983333,2.336566667,
  2.3468,2.357016667,2.367283333,2.377183333,2.387083333,2.39735,2.40725,2.4175,
  2.4274,2.4373,2.44755,2.45745,2.4677,2.477966667,2.488216667,2.498516667,2.5084,
  2.5183,2.5282,2.538466667,2.548366667,2.558266667,2.568516667,2.578783333,
  2.588683333,2.59895,2.6092,2.619466667,2.630066667,2.64065,2.65125,2.662116667,
  2.672716667,2.6833,2.6939,2.7045,2.715083333,2.725683333,2.736266667,2.746866667,
  2.757833333,2.768416667,2.779016667,2.789616667,2.8002,2.810116667,2.820033333,
  2.830316667,2.840216667,2.849766667,2.859316667,2.868866667,2.878783333,2.888683333,
  2.898233333,2.907783333,2.916033333,2.924266667,2.93215,2.940383333,2.947933333,
  2.955816667,2.964066667,2.97195,2.979833333,2.987716667,2.995616667,3.003516667,
  3.011416667,3.01895,3.026833333,3.034366667,3.042266667,3.0498,3.05735,3.065233333,
  3.073133333,3.080666667,3.0882,3.095733333,3.103633333,3.111533333,3.119066667,
  3.126966667,3.134866667,3.14275,3.15065,3.15855,3.166433333,3.174333333,3.182233333,
  3.190133333,3.198016667,3.205916667,3.213166667
};

const vector<double> intensity = {
  1447,2139,1699,755,1258,1070,944,1258,1573,1636,
  1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
  1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
  51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
  184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
  6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
  2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
  2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
};

const vector<double> saturated_pos_min = {
  2.46444988, 2.4746666, 2.4849, 2.49511671, 2.50533342, 2.51556659, 2.52546668, 2.53568339, 2.54563332, 2.55553341, 2.56541657, 2.57566667, 2.58626676, 2.59686661, 2.60778332, 2.61871672, 2.62963343, 2.64056659, 2.6514833, 2.66278338, 2.67406678, 2.68501663, 2.69596672, 2.70693326, 2.71788335, 2.72848344, 2.73943329, 2.75003338, 2.76063323, 2.77121663, 2.78181672, 2.79241657, 2.80299997, 2.8129499, 2.82321668, 2.83313322, 2.84303331, 2.8526001, 2.86213326, 2.87168336, 2.88123322, 2.89078331, 2.9003334, 2.90960002, 2.91886663, 2.92775011, 2.93665004, 2.94589996, 2.95514989, 2.96440005, 2.97364998, 2.9828999, 2.99215007, 3.00139999, 3.01064992, 3.01990008, 3.02915001, 3.03806663, 3.04698324, 3.05591655, 3.0648334, 3.0737834, 3.08270001, 3.09163332, 3.10054994, 3.10946655, 3.1184001, 3.12731671, 3.13623333, 3.14516664, 3.15408325, 3.16300011, 3.17193341, 3.18085003, 3.18976665, 3.19869995, 3.20763326, 3.21623325, 3.22483325, 3.23341656, 3.24201655, 3.25061655, 3.25921679
};

const vector<double> saturated_pos_sec = {
  147.8669928, 148.479996, 149.094, 149.7070026, 150.3200052, 150.93399540000001, 151.5280008, 152.1410034, 152.7379992, 153.33200459999998, 153.92499420000001, 154.5400002, 155.1760056, 155.81199660000001, 156.46699919999998, 157.1230032, 157.77800580000002, 158.43399540000001, 159.088998, 159.7670028, 160.4440068, 161.1009978, 161.7580032, 162.4159956, 163.073001, 163.70900640000002, 164.3659974, 165.0020028, 165.6379938, 166.2729978, 166.9090032, 167.5449942, 168.1799982, 168.776994, 169.39300079999998, 169.9879932, 170.5819986, 171.156006, 171.72799559999999, 172.3010016, 172.8739932, 173.4469986, 174.020004, 174.5760012, 175.13199780000002, 175.6650066, 176.19900239999998, 176.75399760000002, 177.3089934, 177.864003, 178.4189988, 178.973994, 179.5290042, 180.08399939999998, 180.6389952, 181.19400480000002, 181.7490006, 182.28399779999998, 182.81899439999998, 183.354993, 183.890004, 184.427004, 184.9620006, 185.4979992, 186.0329964, 186.567993, 187.104006, 187.6390026, 188.1739998, 188.7099984, 189.244995, 189.7800066, 190.31600459999999, 190.8510018, 191.385999, 191.921997, 192.4579956, 192.973995, 193.489995, 194.0049936, 194.520993, 195.036993, 195.5530074
};

const vector<double> saturated_int = {
  3667.91333, 3829.03906, 3992.62622, 4164.69531, 4438.9165, 4958.67188, 5914.42041, 7855.03125, 11941.041, 21250.4023, 42803.6133, 94525.1094, 216015.453, 472692.219, 961669, 1718756.12, 2641781.25, 3480271.25, 3979093.25, 4087263, 3988863.5, 3942767, 4051667.25, 4250679.5, 4385092, 4301191.5, 3926528, 3335860.5, 2652440.75, 2002597, 1457891.5, 1041989.56, 746353.938, 555698.812, 426332.062, 344843.938, 291420.156, 256906.516, 232813.5, 215017.938, 200963.688, 188910.703, 177631.375, 166692.906, 155609.891, 144831.812, 133511, 121007.531, 108372.43, 97343.3359, 89953.1406, 85699.6328, 82895.1094, 80079.7188, 76910.4375, 73768.9609, 70963.6641, 68590.4766, 66312.75, 64046.4219, 61798.5039, 59813.6211, 58235.5156, 56946.2266, 55774.4766, 54536.9844, 53272.5625, 52113.4883, 51119.7734, 50214.9961, 49111.543, 47641.0234, 45777.4805, 43666.3633, 41646.0508, 39955.7422, 38698.9688, 37738.0273, 36792.582, 35764.0156, 34750.8164, 33958.5547, 33610.1445
};

const vector<double> saturated_cutoff_pos_min = {
  14.3310337, 14.3429499, 14.3551168, 14.3672667, 14.3796835, 14.3923168, 14.4049501, 14.4175835, 14.4299498, 14.4420834, 14.454217, 14.4663496, 14.4782495, 14.4903831, 14.5025167, 14.515317, 14.5275335, 14.5397329, 14.5516834, 14.5636501, 14.5756168, 14.5873337, 14.5993004, 14.6110163, 14.6222496, 14.6334667, 14.6442165, 14.6552162, 14.6661997, 14.6772003, 14.6881838, 14.6987, 14.7094669, 14.7199831, 14.7302504, 14.7405329, 14.7505665, 14.7605829, 14.7703829, 14.7796669, 14.7891998, 14.7985001, 14.810483, 14.8224335, 14.8338833, 14.8455667, 14.8572502, 14.8689499, 14.8806334, 14.8923168, 14.9042501, 14.9164, 14.9287834, 14.9411669, 14.9535503, 14.9659166, 14.9783001, 14.9906836, 15.0025835, 15.0167665, 15.0309162, 15.0450668, 15.0592003, 15.07335, 15.0874996, 15.1016502
};

const vector<double> saturated_cutoff_pos_sec = {
  859.862022, 860.576994, 861.307008, 862.036002, 862.78101, 863.539008, 864.297006, 865.0550099999999, 865.7969879999999, 866.525004, 867.25302, 867.9809759999999, 868.69497, 869.422986, 870.151002, 870.9190199999999, 871.65201, 872.3839740000001, 873.101004, 873.8190060000001, 874.537008, 875.2400220000001, 875.958024, 876.660978, 877.334976, 878.008002, 878.65299, 879.312972, 879.971982, 880.632018, 881.291028, 881.922, 882.5680140000001, 883.198986, 883.815024, 884.431974, 885.03399, 885.6349739999999, 886.222974, 886.780014, 887.351988, 887.910006, 888.62898, 889.3460100000001, 890.032998, 890.7340019999999, 891.4350119999999, 892.1369940000001, 892.8380040000001, 893.539008, 894.255006, 894.9839999999999, 895.7270040000001, 896.470014, 897.213018, 897.9549959999999, 898.698006, 899.441016, 900.15501, 901.00599, 901.854972, 902.704008, 903.5520180000001, 904.401, 905.249976, 906.099012
};

const vector<double> saturated_cutoff_int = {
  1808499.25, 3368120.75, 3803323.25, 4059358, 4092095, 4058075.25, 4160395, 4395341.5, 4573185, 4565417.5, 4371225, 4065336.75, 3716669.5, 3338531, 2982165.25, 2675112, 2464135, 2307864, 2178321, 2053773.88, 1941723.12, 1851044.12, 1768557.75, 1694074, 1643185, 1621222.5, 1615733.12, 1586897.38, 1506051.75, 1374752.25, 1228128.12, 1104698.5, 1004596.56, 922247.438, 846399.25, 775016.375, 717635, 678566.562, 656717.375, 645997, 640481.062, 637876.188, 632300.688, 619110.375, 599831.625, 578701.5, 559403.312, 543297.625, 528049.25, 509227.594, 485023.156, 462585.688, 448690.75, 438401.25, 425110.969, 409034.406, 394403.406, 381114.688, 369241.094, 356554.969, 345889.281, 336646.844, 325948.938, 316341.688, 310746.156, 309641.875
};

const vector<double> cutoff_pos_min = {
  15.34253311, 15.35624981, 15.36995029, 15.38366699, 15.39736652, 15.41156673, 15.42574978, 15.44018364, 15.45436668, 15.46856689, 15.48274994, 15.49695015
};

const vector<double> cutoff_pos_sec = {
  920.5519866, 921.3749885999999, 922.1970174, 923.0200194, 923.8419912, 924.6940038, 925.5449868000001, 926.4110184000001, 927.2620008, 928.1140134, 928.9649964, 929.817009
};

const vector<double> cutoff_int = {
  3.48297429, 15.54384613, 50.31319046, 151.8971405, 411.25631714, 946.44311523, 1642.56152344, 2118.89526367, 2055.13647461, 1665.13232422, 1275.53015137, 1009.70056152
};

MSChromatogram chromatogram;
MSSpectrum spectrum;
for (Size i = 0; i < position.size(); ++i)
{
  chromatogram.push_back(ChromatogramPeak(position[i], intensity[i]));
  spectrum.push_back(Peak1D(position[i], intensity[i]));
}

MSChromatogram saturated_chrom_min, saturated_cutoff_chrom_min, cutoff_chrom_min;
MSChromatogram saturated_chrom_sec, saturated_cutoff_chrom_sec, cutoff_chrom_sec;
MSSpectrum saturated_spec_min, saturated_cutoff_spec_min, cutoff_spec_min;
MSSpectrum saturated_spec_sec, saturated_cutoff_spec_sec, cutoff_spec_sec;
for (Size i = 0; i < saturated_pos_min.size(); ++i)
{
  saturated_chrom_min.push_back(ChromatogramPeak(saturated_pos_min[i], saturated_int[i]));
  saturated_chrom_sec.push_back(ChromatogramPeak(saturated_pos_sec[i], saturated_int[i]));
  saturated_spec_min.push_back(Peak1D(saturated_pos_min[i], saturated_int[i]));
  saturated_spec_sec.push_back(Peak1D(saturated_pos_sec[i], saturated_int[i]));
}
for (Size i = 0; i < saturated_cutoff_pos_min.size(); ++i)
{
  saturated_cutoff_chrom_min.push_back(ChromatogramPeak(saturated_cutoff_pos_min[i], saturated_cutoff_int[i]));
  saturated_cutoff_chrom_sec.push_back(ChromatogramPeak(saturated_cutoff_pos_sec[i], saturated_cutoff_int[i]));
  saturated_cutoff_spec_min.push_back(Peak1D(saturated_cutoff_pos_min[i], saturated_cutoff_int[i]));
  saturated_cutoff_spec_sec.push_back(Peak1D(saturated_cutoff_pos_sec[i], saturated_cutoff_int[i]));
}
for (Size i = 0; i < cutoff_pos_min.size(); ++i)
{
  cutoff_chrom_min.push_back(ChromatogramPeak(cutoff_pos_min[i], cutoff_int[i]));
  cutoff_chrom_sec.push_back(ChromatogramPeak(cutoff_pos_sec[i], cutoff_int[i]));
  cutoff_spec_min.push_back(Peak1D(cutoff_pos_min[i], cutoff_int[i]));
  cutoff_spec_sec.push_back(Peak1D(cutoff_pos_sec[i], cutoff_int[i]));
}


START_SECTION(EmgGradientDescent())
{
  ptr = new EmgGradientDescent();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~EmgGradientDescent())
{
  delete ptr;
}
END_SECTION

START_SECTION(getParameters())
{
  EmgGradientDescent emg;
  const Param& params = emg.getParameters();
  TEST_EQUAL(params.getValue("print_debug"), 0)
  TEST_EQUAL(params.getValue("max_gd_iter"), 100000)
  TEST_EQUAL(params.getValue("compute_additional_points"), "true")
}
END_SECTION

START_SECTION(void fitEMGPeakModel(
  const MSChromatogram& input_peak,
  MSChromatogram& output_peak
) const)
{
  MSChromatogram out_min, out_sec;
  const MSChromatogram::FloatDataArray * fda_emg;

  EmgGradientDescent emg;

  emg.fitEMGPeakModel(chromatogram, out_min);
  TEST_EQUAL(out_min.size(), 107)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 1317410)
  TEST_REAL_SIMILAR((*fda_emg)[1], 2.68121)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0212625)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.0235329)
  // geogebra_print_execute((*fda_emg)[0], (*fda_emg)[1], (*fda_emg)[2], (*fda_emg)[3]);

  emg.fitEMGPeakModel(saturated_chrom_min, out_min);
  emg.fitEMGPeakModel(saturated_chrom_sec, out_sec);
  TEST_EQUAL(out_min.size(), 87)
  TEST_EQUAL(out_min.size(), out_sec.size())
  TOLERANCE_RELATIVE(1.0 + 1e-2)
  for (Size i = 0; i < out_min.size(); i += 9)
  {
    TEST_REAL_SIMILAR(out_min[i].getIntensity(), out_sec[i].getIntensity())
  }
  TOLERANCE_RELATIVE(1.0 + 1e-5)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 7735860)
  TEST_REAL_SIMILAR((*fda_emg)[1], 2.66296)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0394313)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.0394313)
  // geogebra_print_execute((*fda_emg)[0], (*fda_emg)[1], (*fda_emg)[2], (*fda_emg)[3]);

  fda_emg = &out_sec.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 7736020)
  TEST_REAL_SIMILAR((*fda_emg)[1], 159.778)
  TEST_REAL_SIMILAR((*fda_emg)[2], 2.36584)
  TEST_REAL_SIMILAR((*fda_emg)[3], 2.36584)

  emg.fitEMGPeakModel(saturated_cutoff_chrom_min, out_min);
  emg.fitEMGPeakModel(saturated_cutoff_chrom_sec, out_sec);
  TEST_EQUAL(out_min.size(), 71)
  TEST_EQUAL(out_min.size(), out_sec.size())
  TOLERANCE_RELATIVE(1.0 + 1e-2)
  for (Size i = 0; i < out_min.size(); i += 9)
  {
    TEST_REAL_SIMILAR(out_min[i].getIntensity(), out_sec[i].getIntensity())
  }
  TOLERANCE_RELATIVE(1.0 + 1e-5)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 15515900)
  TEST_REAL_SIMILAR((*fda_emg)[1], 14.3453)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0344277)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.188507)
  // geogebra_print_execute((*fda_emg)[0], (*fda_emg)[1], (*fda_emg)[2], (*fda_emg)[3]);

  fda_emg = &out_sec.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 15515900)
  TEST_REAL_SIMILAR((*fda_emg)[1], 860.719)
  TEST_REAL_SIMILAR((*fda_emg)[2], 2.06566)
  TEST_REAL_SIMILAR((*fda_emg)[3], 11.3104)

  emg.fitEMGPeakModel(cutoff_chrom_min, out_min);
  emg.fitEMGPeakModel(cutoff_chrom_sec, out_sec);
  TEST_EQUAL(out_min.size(), 28)
  TEST_EQUAL(out_min.size(), out_sec.size())
  TOLERANCE_RELATIVE(1.0 + 1e-2)
  for (Size i = 0; i < out_min.size(); i += 9)
  {
    TEST_REAL_SIMILAR(out_min[i].getIntensity(), out_sec[i].getIntensity())
  }
  TOLERANCE_RELATIVE(1.0 + 1e-5)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 3791.07)
  TEST_REAL_SIMILAR((*fda_emg)[1], 15.4227)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0210588)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.0476741)
  // geogebra_print_execute((*fda_emg)[0], (*fda_emg)[1], (*fda_emg)[2], (*fda_emg)[3]);

  fda_emg = &out_sec.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 3791.13)
  TEST_REAL_SIMILAR((*fda_emg)[1], 925.363)
  TEST_REAL_SIMILAR((*fda_emg)[2], 1.26351)
  TEST_REAL_SIMILAR((*fda_emg)[3], 2.8605)
}
END_SECTION

START_SECTION(void fitEMGPeakModel(
  const MSSpectrum& input_peak,
  MSSpectrum& output_peak
) const)
{
  MSSpectrum out_min, out_sec;
  const MSSpectrum::FloatDataArray * fda_emg;

  EmgGradientDescent emg;

  emg.fitEMGPeakModel(spectrum, out_min);
  TEST_EQUAL(out_min.size(), 107)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 1317410)
  TEST_REAL_SIMILAR((*fda_emg)[1], 2.68121)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0212625)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.0235329)

  emg.fitEMGPeakModel(saturated_spec_min, out_min);
  emg.fitEMGPeakModel(saturated_spec_sec, out_sec);
  TEST_EQUAL(out_min.size(), 87)
  TEST_EQUAL(out_min.size(), out_sec.size())
  TOLERANCE_RELATIVE(1.0 + 1e-2)
  for (Size i = 0; i < out_min.size(); i += 9)
  {
    TEST_REAL_SIMILAR(out_min[i].getIntensity(), out_sec[i].getIntensity())
  }
  TOLERANCE_RELATIVE(1.0 + 1e-5)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 7735860)
  TEST_REAL_SIMILAR((*fda_emg)[1], 2.66296)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0394313)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.0394313)

  fda_emg = &out_sec.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 7736020)
  TEST_REAL_SIMILAR((*fda_emg)[1], 159.778)
  TEST_REAL_SIMILAR((*fda_emg)[2], 2.36584)
  TEST_REAL_SIMILAR((*fda_emg)[3], 2.36584)

  emg.fitEMGPeakModel(saturated_cutoff_spec_min, out_min);
  emg.fitEMGPeakModel(saturated_cutoff_spec_sec, out_sec);
  TEST_EQUAL(out_min.size(), 71)
  TEST_EQUAL(out_min.size(), out_sec.size())
  TOLERANCE_RELATIVE(1.0 + 1e-2)
  for (Size i = 0; i < out_min.size(); i += 9)
  {
    TEST_REAL_SIMILAR(out_min[i].getIntensity(), out_sec[i].getIntensity())
  }
  TOLERANCE_RELATIVE(1.0 + 1e-5)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 15515900)
  TEST_REAL_SIMILAR((*fda_emg)[1], 14.3453)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0344277)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.188507)

  fda_emg = &out_sec.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 15515900)
  TEST_REAL_SIMILAR((*fda_emg)[1], 860.719)
  TEST_REAL_SIMILAR((*fda_emg)[2], 2.06566)
  TEST_REAL_SIMILAR((*fda_emg)[3], 11.3104)

  emg.fitEMGPeakModel(cutoff_spec_min, out_min);
  emg.fitEMGPeakModel(cutoff_spec_sec, out_sec);
  TEST_EQUAL(out_min.size(), 28)
  TEST_EQUAL(out_min.size(), out_sec.size())
  TOLERANCE_RELATIVE(1.0 + 1e-2)
  for (Size i = 0; i < out_min.size(); i += 9)
  {
    TEST_REAL_SIMILAR(out_min[i].getIntensity(), out_sec[i].getIntensity())
  }
  TOLERANCE_RELATIVE(1.0 + 1e-5)

  fda_emg = &out_min.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 3791.07)
  TEST_REAL_SIMILAR((*fda_emg)[1], 15.4227)
  TEST_REAL_SIMILAR((*fda_emg)[2], 0.0210588)
  TEST_REAL_SIMILAR((*fda_emg)[3], 0.0476741)

  fda_emg = &out_sec.getFloatDataArrays()[0];
  TEST_EQUAL(fda_emg->getName(), "emg_parameters")
  TEST_REAL_SIMILAR((*fda_emg)[0], 3791.13)
  TEST_REAL_SIMILAR((*fda_emg)[1], 925.363)
  TEST_REAL_SIMILAR((*fda_emg)[2], 1.26351)
  TEST_REAL_SIMILAR((*fda_emg)[3], 2.8605)
}
END_SECTION

START_SECTION(double Loss_function(
  const std::vector<double>& xs,
  const std::vector<double>& ys,
  const double h,
  const double mu,
  const double sigma,
  const double tau
) const)
{
  MSChromatogram out_min;
  EmgGradientDescent emg;
  Param params = emg.getParameters();
  params.setValue("compute_additional_points", "false");
  emg.setParameters(params);
  emg.fitEMGPeakModel(chromatogram, out_min);
  EmgGradientDescent_friend emg_f;

  const MSChromatogram::FloatDataArray& fda_emg = out_min.getFloatDataArrays()[0];
  TEST_REAL_SIMILAR(emg_f.Loss_function(position, intensity, fda_emg[0], fda_emg[1], fda_emg[2], fda_emg[3]), 60778399.8312241)
  // geogebra_print_execute(fda_emg[0], fda_emg[1], fda_emg[2], fda_emg[3]);

  emg.fitEMGPeakModel(saturated_chrom_min, out_min);
  const MSChromatogram::FloatDataArray& fda_emg_sat = out_min.getFloatDataArrays()[0];
  TEST_REAL_SIMILAR(emg_f.Loss_function(saturated_pos_min, saturated_int, fda_emg_sat[0], fda_emg_sat[1], fda_emg_sat[2], fda_emg_sat[3]), 187412764882.422)
  // geogebra_print_execute(fda_emg_sat[0], fda_emg_sat[1], fda_emg_sat[2], fda_emg_sat[3]);

  emg.fitEMGPeakModel(saturated_cutoff_chrom_min, out_min);
  const MSChromatogram::FloatDataArray& fda_emg_sat_cut = out_min.getFloatDataArrays()[0];
  TEST_REAL_SIMILAR(emg_f.Loss_function(saturated_cutoff_pos_min, saturated_cutoff_int, fda_emg_sat_cut[0], fda_emg_sat_cut[1], fda_emg_sat_cut[2], fda_emg_sat_cut[3]), 56213636966.189)
  // geogebra_print_execute(fda_emg_sat_cut[0], fda_emg_sat_cut[1], fda_emg_sat_cut[2], fda_emg_sat_cut[3]);

  emg.fitEMGPeakModel(cutoff_chrom_min, out_min);
  const MSChromatogram::FloatDataArray& fda_emg_cut = out_min.getFloatDataArrays()[0];
  TEST_REAL_SIMILAR(emg_f.Loss_function(cutoff_pos_min, cutoff_int, fda_emg_cut[0], fda_emg_cut[1], fda_emg_cut[2], fda_emg_cut[3]), 651.824632922326)
  // geogebra_print_execute(fda_emg_cut[0], fda_emg_cut[1], fda_emg_cut[2], fda_emg_cut[3]);
}
END_SECTION

START_SECTION(void extractTrainingSet(
  const std::vector<double>& xs,
  const std::vector<double>& ys,
  std::vector<double>& TrX,
  std::vector<double>& TrY
) const)
{
  EmgGradientDescent_friend emg_f;
  vector<double> TrX, TrY;

  emg_f.extractTrainingSet(position, intensity, TrX, TrY);
  // non-saturated, non-cutoff peak: no point is filtered, all are valid
  TEST_EQUAL(TrX.size(), position.size())

  emg_f.extractTrainingSet(saturated_pos_min, saturated_int, TrX, TrY);
  TEST_NOT_EQUAL(TrX.size(), saturated_pos_min.size())
  TEST_EQUAL(TrX.size(), 77)

  emg_f.extractTrainingSet(saturated_cutoff_pos_min, saturated_cutoff_int, TrX, TrY);
  TEST_NOT_EQUAL(TrX.size(), saturated_cutoff_pos_min.size())
  TEST_EQUAL(TrX.size(), 61)

  emg_f.extractTrainingSet(cutoff_pos_min, cutoff_int, TrX, TrY);
  // cutoff but non-saturated peak: no point is filtered, all are valid
  TEST_EQUAL(TrX.size(), cutoff_pos_min.size())
}
END_SECTION

START_SECTION(double computeMuMaxDistance(const std::vector<double>& xs) const)
{
  EmgGradientDescent_friend emg_f;
  vector<double> xs { 3, 2, 4, 2, 4, 5, 7, 9, 3 };
  TEST_REAL_SIMILAR(emg_f.computeMuMaxDistance(xs), 2.45)
  xs.clear();
  TEST_REAL_SIMILAR(emg_f.computeMuMaxDistance(xs), 0.0) // empty vector case
}
END_SECTION

START_SECTION(double computeInitialMean(
  const std::vector<double>& xs,
  const std::vector<double>& ys
) const)
{
  EmgGradientDescent_friend emg_f;
  double mu;
  mu = emg_f.computeInitialMean(position, intensity);
  TEST_REAL_SIMILAR(mu, 2.69743333333333)
  mu = emg_f.computeInitialMean(saturated_pos_min, saturated_int);
  TEST_REAL_SIMILAR(mu, 2.69516110583333)
  mu = emg_f.computeInitialMean(saturated_cutoff_pos_sec, saturated_cutoff_int);
  TEST_REAL_SIMILAR(mu, 865.1314205)
  mu = emg_f.computeInitialMean(cutoff_pos_sec, cutoff_int);
  TEST_REAL_SIMILAR(mu, 926.90050115)
  const vector<double> empty;
  TEST_EXCEPTION(Exception::SizeUnderflow, emg_f.computeInitialMean(empty, empty))
}
END_SECTION

START_SECTION(void iRpropPlus(
  const double prev_diff_E_param,
  double& diff_E_param,
  double& param_lr,
  double& param_update,
  double& param,
  const double current_E,
  const double previous_E
) const)
{
  EmgGradientDescent_friend emg_f;
  const double prev_diff_E_param { 10.0 };
  double diff_E_param { 20.0 };
  double param_lr { 4.0 };
  double param_update { 0.5 };
  double param { 860.0 };
  const double current_E { 13.0 };
  const double previous_E { 14.0 };

  emg_f.iRpropPlus(
    prev_diff_E_param, diff_E_param, param_lr,
    param_update, param, current_E, previous_E
  );
  TEST_REAL_SIMILAR(diff_E_param, 20.0)
  TEST_REAL_SIMILAR(param_lr, 4.8)
  TEST_REAL_SIMILAR(param_update, -4.8)
  TEST_REAL_SIMILAR(param, 855.2)

  diff_E_param = -20.0;
  param_lr = 4.0;
  param_update = 0.5;
  param = 860.0;
  emg_f.iRpropPlus(
    prev_diff_E_param, diff_E_param, param_lr,
    param_update, param, current_E, previous_E
  );
  TEST_REAL_SIMILAR(diff_E_param, 0.0)
  TEST_REAL_SIMILAR(param_lr, 2.0)
  TEST_REAL_SIMILAR(param_update, 0.5)
  TEST_REAL_SIMILAR(param, 860.0)

  diff_E_param = 0.0;
  param_lr = 4.0;
  param_update = 0.5;
  param = 860.0;
  emg_f.iRpropPlus(
    prev_diff_E_param, diff_E_param, param_lr,
    param_update, param, current_E, previous_E
  );
  TEST_REAL_SIMILAR(diff_E_param, 0.0)
  TEST_REAL_SIMILAR(param_lr, 4.0)
  TEST_REAL_SIMILAR(param_update, -4.0)
  TEST_REAL_SIMILAR(param, 856.0)
}
END_SECTION

START_SECTION(double compute_z(
  const double x,
  const double mu,
  const double sigma,
  const double tau
) const)
{
  EmgGradientDescent_friend emg_f;
  double x;
  const double mu { 14.3453 };
  const double sigma { 0.0344277 };
  const double tau { 0.188507 };

  x = mu - 1.0 / 60.0;
  TEST_REAL_SIMILAR(emg_f.compute_z(x, mu, sigma, tau), 0.471456263584609)

  x = mu + 1.0 / 60.0;
  TEST_REAL_SIMILAR(emg_f.compute_z(x, mu, sigma, tau), -0.213173439809831)

  x = -3333333;
  TEST_REAL_SIMILAR(emg_f.compute_z(x, mu, sigma, tau), 68463258.2588395)
}
END_SECTION

START_SECTION(double emg_point(
  const double x,
  const double h,
  const double mu,
  const double sigma,
  const double tau
) const)
{
  EmgGradientDescent_friend emg_f;
  double x;
  double h { 15515900 };
  double mu { 14.3453 };
  double sigma { 0.0344277 };
  double tau { 0.188507 };

  x = mu - 1.0 / 60.0;
  TEST_REAL_SIMILAR(emg_f.emg_point(x, h, mu, sigma, tau), 1992032.65711041)

  x = mu + 1.0 / 60.0;
  TEST_REAL_SIMILAR(emg_f.emg_point(x, h, mu, sigma, tau), 4088964.97520213)

  x = -3333333;
  TEST_REAL_SIMILAR(emg_f.emg_point(x, h, mu, sigma, tau), 0.0)

  mu = 860.719;
  sigma = 2.06566;
  tau = 11.3104;

  x = mu - 1;
  TEST_REAL_SIMILAR(emg_f.emg_point(x, h, mu, sigma, tau), 1992033.06584247)

  x = mu + 1;
  TEST_REAL_SIMILAR(emg_f.emg_point(x, h, mu, sigma, tau), 4088968.52957875)

  x = -200000000;
  TEST_REAL_SIMILAR(emg_f.emg_point(x, h, mu, sigma, tau), 0.0)
}
END_SECTION

START_SECTION(void applyEstimatedParameters(
  const std::vector<double>& xs,
  const double h,
  const double mu,
  const double sigma,
  const double tau,
  std::vector<double>& out_xs,
  std::vector<double>& out_ys
) const)
{
  EmgGradientDescent_friend emg_f;
  const double h { 15515900 };
  const double mu { 14.3453 };
  const double sigma { 0.0344277 };
  const double tau { 0.188507 };
  vector<double> out_xs;
  vector<double> out_ys;

  Param params = emg_f.emg_gd_.getParameters();

  params.setValue("compute_additional_points", "false");
  emg_f.emg_gd_.setParameters(params);
  emg_f.applyEstimatedParameters(saturated_cutoff_pos_min, h, mu, sigma, tau, out_xs, out_ys);
  TEST_EQUAL(out_xs.size(), saturated_cutoff_pos_min.size())
  TEST_REAL_SIMILAR(out_xs.front(), 14.3310337)
  TEST_REAL_SIMILAR(out_ys.front(), 2144281.1472228)

  params.setValue("compute_additional_points", "true");
  emg_f.emg_gd_.setParameters(params);
  emg_f.applyEstimatedParameters(saturated_cutoff_pos_min, h, mu, sigma, tau, out_xs, out_ys);
  TEST_EQUAL(out_xs.size(), 71) // more points than before
  TEST_REAL_SIMILAR(out_xs.front(), 14.2717555076923) // peak was cutoff on the left side
  TEST_REAL_SIMILAR(out_ys.front(), 108845.941990663)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
