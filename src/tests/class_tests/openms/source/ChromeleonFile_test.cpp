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
#include <OpenMS/FORMAT/ChromeleonFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/File.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ChromeleonFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromeleonFile* ptr = 0;
ChromeleonFile* null_ptr = 0;

START_SECTION(ChromeleonFile())
{
  ptr = new ChromeleonFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ChromeleonFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& filename, MSExperiment& experiment) const)
{
  String input_filepath = OPENMS_GET_TEST_DATA_PATH("20171013_HMP_C61_ISO_P1_GA1_UV_VIS_2.txt");
  MSExperiment experiment;
  ChromeleonFile cf;
  cf.load(input_filepath, experiment);
  TEST_EQUAL(experiment.getMetaValue("acq_method_name"), "UV_VIS_2")
  TEST_EQUAL(experiment.getMetaValue("mzml_id"), "20171013_C61_ISO_P1_GA1")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getName(), "HM_metode_ZorBax_0,02%_Acetic_acid_ver6")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getSoftware().getName(), "New ProcMethod")
  TEST_EQUAL(experiment.getMetaValue("injection_date"), "10/13/2017")
  TEST_EQUAL(experiment.getMetaValue("injection_time"), "6:28:26 PM")
  TEST_EQUAL(experiment.getMetaValue("detector"), "UV")
  TEST_EQUAL(experiment.getMetaValue("signal_quantity"), "Absorbance")
  TEST_EQUAL(experiment.getMetaValue("signal_unit"), "mAU")
  TEST_EQUAL(experiment.getMetaValue("signal_info"), "WVL:280 nm")
  const vector<MSChromatogram> chromatograms = experiment.getChromatograms();
  TEST_EQUAL(chromatograms.size(), 1);
  TEST_EQUAL(chromatograms[0].size(), 3301);
  const MSChromatogram& c = chromatograms[0];
  TEST_REAL_SIMILAR(c[0].getRT(), 0.0)
  TEST_REAL_SIMILAR(c[0].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(c[660].getRT(), 2.2)
  TEST_REAL_SIMILAR(c[660].getIntensity(), -0.812998)
  TEST_REAL_SIMILAR(c[1320].getRT(), 4.4)
  TEST_REAL_SIMILAR(c[1320].getIntensity(), -0.791189)
  TEST_REAL_SIMILAR(c[1980].getRT(), 6.6)
  TEST_REAL_SIMILAR(c[1980].getIntensity(), -0.285533)
  TEST_REAL_SIMILAR(c[2640].getRT(), 8.8)
  TEST_REAL_SIMILAR(c[2640].getIntensity(), -0.485941)
  TEST_REAL_SIMILAR(c[3300].getRT(), 11.0)
  TEST_REAL_SIMILAR(c[3300].getIntensity(), -0.130904)

  MzMLFile mzml;
  const String output_filepath = File::getTemporaryFile();
  mzml.store(output_filepath, experiment);
  MSExperiment read_exp;
  mzml.load(output_filepath, read_exp);
  TEST_EQUAL(read_exp.getChromatograms().size(), 1);
  const MSChromatogram& c1 = experiment.getChromatograms()[0];
  const MSChromatogram& c2 = read_exp.getChromatograms()[0];
  TEST_EQUAL(c1.size(), c2.size())
  TEST_REAL_SIMILAR(c1[0].getRT(), c2[0].getRT())
  TEST_REAL_SIMILAR(c1[0].getIntensity(), c2[0].getIntensity())
  TEST_REAL_SIMILAR(c1[660].getRT(), c2[660].getRT())
  TEST_REAL_SIMILAR(c1[660].getIntensity(), c2[660].getIntensity())
  TEST_REAL_SIMILAR(c1[1320].getRT(), c2[1320].getRT())
  TEST_REAL_SIMILAR(c1[1320].getIntensity(), c2[1320].getIntensity())
  TEST_REAL_SIMILAR(c1[1980].getRT(), c2[1980].getRT())
  TEST_REAL_SIMILAR(c1[1980].getIntensity(), c2[1980].getIntensity())
  TEST_REAL_SIMILAR(c1[2640].getRT(), c2[2640].getRT())
  TEST_REAL_SIMILAR(c1[2640].getIntensity(), c2[2640].getIntensity())
  TEST_REAL_SIMILAR(c1[3300].getRT(), c2[3300].getRT())
  TEST_REAL_SIMILAR(c1[3300].getIntensity(), c2[3300].getIntensity())
}
END_SECTION

START_SECTION(load_with_new_raw_data_header)
{
  String input_filepath = OPENMS_GET_TEST_DATA_PATH("ChromeleonFile_new_header.txt");
  MSExperiment experiment;
  ChromeleonFile cf;
  cf.load(input_filepath, experiment);
  TEST_EQUAL(experiment.getMetaValue("acq_method_name"), "RID_Signal")
  TEST_EQUAL(experiment.getMetaValue("mzml_id"), "S1")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getName(), "SUGARS_MP.M")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getSoftware().getName(), "SUGARS_CAL")
  TEST_EQUAL(experiment.getMetaValue("injection_date"), "13/06/2019")
  TEST_EQUAL(experiment.getMetaValue("injection_time"), "12:11:41 AM")
  TEST_EQUAL(experiment.getMetaValue("detector"), "LCSystem")
  TEST_EQUAL(experiment.getMetaValue("signal_quantity"), "")
  TEST_EQUAL(experiment.getMetaValue("signal_unit"), "nRIU")
  TEST_EQUAL(experiment.getMetaValue("signal_info"), "")
  const vector<MSChromatogram> chromatograms = experiment.getChromatograms();
  TEST_EQUAL(chromatograms.size(), 1);
  TEST_EQUAL(chromatograms[0].size(), 10);
  const MSChromatogram& c = chromatograms[0];
  TEST_REAL_SIMILAR(c[0].getRT(), 0.0)
  TEST_REAL_SIMILAR(c[0].getIntensity(), 5.060000)
  TEST_REAL_SIMILAR(c[2].getRT(), 0.014430)
  TEST_REAL_SIMILAR(c[2].getIntensity(), 5.450000)
  TEST_REAL_SIMILAR(c[4].getRT(), 0.028860)
  TEST_REAL_SIMILAR(c[4].getIntensity(), 5.580000)
  TEST_REAL_SIMILAR(c[6].getRT(), 0.043290)
  TEST_REAL_SIMILAR(c[6].getIntensity(), 5.380000)
  TEST_REAL_SIMILAR(c[9].getRT(), 0.064935)
  TEST_REAL_SIMILAR(c[9].getIntensity(), 4.930000)

  MzMLFile mzml;
  const String output_filepath = File::getTemporaryFile();
  mzml.store(output_filepath, experiment);
  MSExperiment read_exp;
  mzml.load(output_filepath, read_exp);
  TEST_EQUAL(read_exp.getChromatograms().size(), 1);
  const MSChromatogram& c1 = experiment.getChromatograms()[0];
  const MSChromatogram& c2 = read_exp.getChromatograms()[0];
  TEST_EQUAL(c1.size(), c2.size())
  TEST_REAL_SIMILAR(c1[0].getRT(), c2[0].getRT())
  TEST_REAL_SIMILAR(c1[0].getIntensity(), c2[0].getIntensity())
  TEST_REAL_SIMILAR(c1[2].getRT(), c2[2].getRT())
  TEST_REAL_SIMILAR(c1[2].getIntensity(), c2[2].getIntensity())
  TEST_REAL_SIMILAR(c1[4].getRT(), c2[4].getRT())
  TEST_REAL_SIMILAR(c1[4].getIntensity(), c2[4].getIntensity())
  TEST_REAL_SIMILAR(c1[6].getRT(), c2[6].getRT())
  TEST_REAL_SIMILAR(c1[6].getIntensity(), c2[6].getIntensity())
  TEST_REAL_SIMILAR(c1[9].getRT(), c2[9].getRT())
  TEST_REAL_SIMILAR(c1[9].getIntensity(), c2[9].getIntensity())
}
END_SECTION

START_SECTION(load_file_with_comma_thousands_separator)
{
  String input_filepath = OPENMS_GET_TEST_DATA_PATH("ChromeleonFile_commas.txt");
  MSExperiment experiment;
  ChromeleonFile cf;
  cf.load(input_filepath, experiment);
  TEST_EQUAL(experiment.getMetaValue("acq_method_name"), "RID_Signal")
  TEST_EQUAL(experiment.getMetaValue("mzml_id"), "S2")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getName(), "SUGARS_MP.M")
  TEST_EQUAL(experiment.getExperimentalSettings().getInstrument().getSoftware().getName(), "SUGARS_CAL")
  TEST_EQUAL(experiment.getMetaValue("injection_date"), "12/06/2019")
  TEST_EQUAL(experiment.getMetaValue("injection_time"), "11:49:36 PM")
  TEST_EQUAL(experiment.getMetaValue("detector"), "LCSystem")
  TEST_EQUAL(experiment.getMetaValue("signal_quantity"), "")
  TEST_EQUAL(experiment.getMetaValue("signal_unit"), "nRIU")
  TEST_EQUAL(experiment.getMetaValue("signal_info"), "")
  const vector<MSChromatogram> chromatograms = experiment.getChromatograms();
  TEST_EQUAL(chromatograms.size(), 1);
  TEST_EQUAL(chromatograms[0].size(), 8);
  const MSChromatogram& c = chromatograms[0];
  TEST_REAL_SIMILAR(c[0].getRT(), 0.0)
  TEST_REAL_SIMILAR(c[1].getRT(), 8.300000)
  TEST_REAL_SIMILAR(c[2].getRT(), 8.513709)
  TEST_REAL_SIMILAR(c[3].getRT(), 8.520924)
  TEST_REAL_SIMILAR(c[4].getRT(), 9.855700)
  TEST_REAL_SIMILAR(c[5].getRT(), 9.884560)
  TEST_REAL_SIMILAR(c[6].getRT(), 9.898991)
  TEST_REAL_SIMILAR(c[7].getRT(), 9.920635)
  TEST_REAL_SIMILAR(c[0].getIntensity(), 1.4)
  TEST_REAL_SIMILAR(c[1].getIntensity(), 1.6)
  TEST_REAL_SIMILAR(c[2].getIntensity(), -18.980000)
  TEST_REAL_SIMILAR(c[3].getIntensity(), -1234567.890000)
  TEST_REAL_SIMILAR(c[4].getIntensity(), 1946.610000)
  TEST_REAL_SIMILAR(c[5].getIntensity(), 2067.450000)
  TEST_REAL_SIMILAR(c[6].getIntensity(), 2345678.900000)
  TEST_REAL_SIMILAR(c[7].getIntensity(), 2028.580000)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
