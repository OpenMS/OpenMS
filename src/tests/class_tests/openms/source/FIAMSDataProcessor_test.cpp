// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FIAMSDataProcessor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FIAMSDataProcessor* ptr = nullptr;
FIAMSDataProcessor* null_ptr = nullptr;
START_SECTION(FIAMSDataProcessor())
{
    ptr = new FIAMSDataProcessor();
    TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(virtual ~FIAMSDataProcessor())
{
    delete ptr;
}
END_SECTION

String filename = "SerumTest";

File::TempDir temp_dir;

FIAMSDataProcessor fia_processor;
Param p;
p.setValue("filename", filename);
p.setValue("dir_output", temp_dir.getPath());
p.setValue("resolution", 120000.0);
p.setValue("polarity", "negative");
p.setValue("max_mz", 1500);
p.setValue("bin_step", 20);
p.setValue("db:mapping", std::vector<std::string>{OPENMS_GET_TEST_DATA_PATH("reducedHMDBMapping.tsv")});
p.setValue("db:struct", std::vector<std::string>{OPENMS_GET_TEST_DATA_PATH("reducedHMDB2StructMapping.tsv")});
p.setValue("positive_adducts", OPENMS_GET_TEST_DATA_PATH("FIAMS_negative_adducts.tsv"));
p.setValue("negative_adducts", OPENMS_GET_TEST_DATA_PATH("FIAMS_positive_adducts.tsv"));
fia_processor.setParameters(p);

MSExperiment exp;
MzMLFile mzml;
mzml.load(String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input")) + "/" + filename + ".mzML", exp);

MSExperiment exp_merged;
MzMLFile mzml_merged;
mzml.load(String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input")) + "/" + filename + "_merged.mzML", exp_merged);
MSSpectrum spec_merged = exp_merged.getSpectra()[0];

MSExperiment exp_picked;
MzMLFile mzml_picked;
mzml.load(String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input")) + "/" + filename + "_picked.mzML", exp_picked);
MSSpectrum spec_picked = exp_picked.getSpectra()[0];

PeakMap input;
Peak1D peak;
std::vector<float> ints {100, 120, 130, 140, 150, 100, 60, 50, 30};
std::vector<float> rts {10, 20, 30, 40};
std::vector<MSSpectrum> spectra;
for (Size i = 0; i < rts.size(); ++i) {
    MSSpectrum s;
    for (Size j = 0; j < ints.size(); ++j) {
        peak.setIntensity(ints[j]); peak.setMZ(100 + j*2);
        s.push_back(peak);
    }
    s.setRT(rts[i]);
    input.addSpectrum(s);
    spectra.push_back(s);
}
MSSpectrum merged = fia_processor.mergeAlongTime(spectra);

START_SECTION((void cutForTime(const MSExperiment & experiment, vector<MSSpectrum> & output, float n_seconds)))
{
    vector<MSSpectrum> output1;
    fia_processor.cutForTime(input, 0, output1);
    TEST_EQUAL(output1.size(), 0);
    vector<MSSpectrum> output2;
    fia_processor.cutForTime(input, 25, output2);
    TEST_EQUAL(output2.size(), 2);
    vector<MSSpectrum> output3;
    fia_processor.cutForTime(input, 100, output3);
    TEST_EQUAL(output3.size(), 4);
    PeakMap empty_input;
    vector<MSSpectrum> output4;
    fia_processor.cutForTime(empty_input, 100, output4);
    TEST_EQUAL(output4.size(), 0);
}
END_SECTION

START_SECTION((mergeAlongTime))
{
    MSSpectrum output = fia_processor.mergeAlongTime(spectra);
    TEST_EQUAL(!output.empty(), true);
    TEST_EQUAL(abs(output.MZBegin(100)->getIntensity() - 400.0) < 1, true);
    TEST_EQUAL(abs(output.MZBegin(102)->getIntensity() - 480.0) < 1, true);
}
END_SECTION

START_SECTION((extractPeaks))
{
    MSSpectrum picked = fia_processor.extractPeaks(merged);
    TEST_EQUAL(abs(picked.MZBegin(108)->getIntensity() - 133) < 1, true);
    TEST_EQUAL(abs(picked.MZBegin(112)->getIntensity() - 66) < 1, true);
}
END_SECTION

START_SECTION((convertToFeatureMap))
{
    MSSpectrum picked;
    for (Size j = 0; j < 10; ++j) {
        peak.setIntensity(50); peak.setMZ(100 + j*2);
        picked.push_back(peak);
    }
    FeatureMap output_feature = fia_processor.convertToFeatureMap(picked);
    for (auto it = output_feature.begin(); it != output_feature.end(); ++it) {
        TEST_EQUAL(it->getIntensity() == 50, true);
    }
}
END_SECTION

START_SECTION((test_run_cached))
{
    MzTab mztab_output_30;
    fia_processor.run(exp, 30, mztab_output_30);
    String filename_30 = "SerumTest_merged_30.mzML";
    TEST_EQUAL(File::exists(temp_dir.getPath() + filename_30), true);
    bool is_cached_after = fia_processor.run(exp, 30, mztab_output_30);
    TEST_EQUAL(is_cached_after, true);
}
END_SECTION

START_SECTION((test_run_empty))
{
    MzTab mztab_output_0;
    String filename_0 = "SerumTest_picked_0.mzML";
    String filename_mztab = "SerumTest_0.mzTab";
    fia_processor.run(exp, 0, mztab_output_0);
    TEST_EQUAL(File::exists(temp_dir.getPath() + filename_0), true);
    TEST_EQUAL(File::exists(temp_dir.getPath() + filename_mztab), true);
    TEST_EQUAL(mztab_output_0.getPSMSectionRows().size(), 0);
}
END_SECTION


START_SECTION((test_run_full))
{
    vector<MSSpectrum> spec_vec;
    fia_processor.cutForTime(exp, 1000, spec_vec);
    MSSpectrum merged_result = fia_processor.mergeAlongTime(spec_vec);
    vector<double> mzs {109.951239, 109.962281, 109.986031, 109.999156};
    for (double mz : mzs) {
        TEST_REAL_SIMILAR(merged_result.MZBegin(mz)->getIntensity(), spec_merged.MZBegin(mz)->getIntensity());
    }
    MSSpectrum picked_result = fia_processor.extractPeaks(merged_result);
    vector<double> mzs_picked {109.951246, 109.957552, 109.959885, 109.961982, 109.982828};
    for (double mz : mzs_picked) {
        TEST_EQUAL(picked_result.MZBegin(mz)->getIntensity(), spec_picked.MZBegin(mz)->getIntensity());
    }
    TEST_EQUAL(picked_result.size(), spec_picked.size());
    MzTab mztab_output;
    fia_processor.run(exp, 1000, mztab_output);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST