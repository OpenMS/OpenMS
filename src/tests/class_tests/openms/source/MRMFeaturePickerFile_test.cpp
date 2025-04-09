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
#include <OpenMS/FORMAT/MRMFeaturePickerFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMFeaturePickerFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeaturePickerFile* ptr = nullptr;
MRMFeaturePickerFile* null_ptr = nullptr;
const String filepath = OPENMS_GET_TEST_DATA_PATH("MRMFeaturePickerFile.csv");

START_SECTION(MRMFeaturePickerFile())
{
  ptr = new MRMFeaturePickerFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MRMFeaturePickerFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(void load(
  const String& filename,
  std::vector<MRMFeaturePicker::ComponentParams>& cp_list,
  std::vector<MRMFeaturePicker::ComponentGroupParams>& cgp_list
))
{
  MRMFeaturePickerFile file;
  std::vector<MRMFeaturePicker::ComponentParams> cp_list;
  std::vector<MRMFeaturePicker::ComponentGroupParams> cgp_list;
  file.load(filepath, cp_list, cgp_list);

  TEST_EQUAL(cp_list.size(), 11)
  TEST_EQUAL(cgp_list.size(), 5)

  TEST_EQUAL(cp_list[1].component_name, "arg-L.arg-L_1.Light")
  TEST_EQUAL(cp_list[1].component_group_name, "arg-L")
  TEST_EQUAL(cp_list[1].params.getValue("sgolay_frame_length"), 152)
  TEST_EQUAL(cp_list[1].params.getValue("sgolay_polynomial_order"), 32)
  TEST_REAL_SIMILAR(cp_list[1].params.getValue("gauss_width"), 0.152)
  TEST_EQUAL(cp_list[1].params.getValue("use_gauss"), "false")
  TEST_REAL_SIMILAR(cp_list[1].params.getValue("peak_width"), 0.12)
  TEST_REAL_SIMILAR(cp_list[1].params.getValue("signal_to_noise"), 0.012)
  TEST_REAL_SIMILAR(cp_list[1].params.getValue("sn_win_len"), 10002.0)
  TEST_EQUAL(cp_list[1].params.getValue("sn_bin_count"), 302)
  TEST_EQUAL(cp_list[1].params.getValue("write_sn_log_messages"), "false")
  TEST_EQUAL(cp_list[1].params.getValue("remove_overlapping_peaks"), "true")
  TEST_EQUAL(cp_list[1].params.getValue("method"), "corrected2")

  TEST_EQUAL(cp_list[9].component_name, "ser-L.ser-L_2.Light")
  TEST_EQUAL(cp_list[9].component_group_name, "ser-L")
  TEST_EQUAL(cp_list[9].params.getValue("sgolay_frame_length"), 160)
  TEST_EQUAL(cp_list[9].params.getValue("sgolay_polynomial_order"), 40)
  TEST_REAL_SIMILAR(cp_list[9].params.getValue("gauss_width"), 0.16)
  TEST_EQUAL(cp_list[9].params.getValue("use_gauss"), "false")
  TEST_REAL_SIMILAR(cp_list[9].params.getValue("peak_width"), 0.2)
  TEST_REAL_SIMILAR(cp_list[9].params.getValue("signal_to_noise"), 0.02)
  TEST_REAL_SIMILAR(cp_list[9].params.getValue("sn_win_len"), 10010.0)
  TEST_EQUAL(cp_list[9].params.getValue("sn_bin_count"), 310)
  TEST_EQUAL(cp_list[9].params.getValue("write_sn_log_messages"), "false")
  TEST_EQUAL(cp_list[9].params.getValue("remove_overlapping_peaks"), "true")
  TEST_EQUAL(cp_list[9].params.getValue("method"), "corrected10")

  TEST_EQUAL(cp_list[10].component_name, "component2")
  TEST_EQUAL(cp_list[10].component_group_name, "group2")

  TEST_EQUAL(cp_list[10].params.getValue("sgolay_polynomial_order"), 43)
  TEST_REAL_SIMILAR(cp_list[10].params.getValue("gauss_width"), 0.163)
  TEST_EQUAL(cp_list[10].params.getValue("use_gauss"), "true")
  TEST_REAL_SIMILAR(cp_list[10].params.getValue("peak_width"), 0.23)
  TEST_REAL_SIMILAR(cp_list[10].params.getValue("signal_to_noise"), 0.023)
  TEST_REAL_SIMILAR(cp_list[10].params.getValue("sn_win_len"), 10013.0)
  TEST_EQUAL(cp_list[10].params.getValue("sn_bin_count"), 313)
  TEST_EQUAL(cp_list[10].params.getValue("write_sn_log_messages"), "true")
  TEST_EQUAL(cp_list[10].params.getValue("remove_overlapping_peaks"), "false")
  TEST_EQUAL(cp_list[10].params.getValue("method"), "corrected13")

  TEST_EQUAL(cp_list[10].params.exists("sgolay_frame_length"), false)

  TEST_EQUAL(cgp_list[1].component_group_name, "orn")
  TEST_EQUAL(cgp_list[1].params.getValue("stop_after_feature"), 6)
  TEST_REAL_SIMILAR(cgp_list[1].params.getValue("stop_after_intensity_ratio"), 0.0006)
  TEST_REAL_SIMILAR(cgp_list[1].params.getValue("min_peak_width"), -6.0)
  TEST_EQUAL(cgp_list[1].params.getValue("peak_integration"), "smoothed3")
  TEST_EQUAL(cgp_list[1].params.getValue("background_subtraction"), "none3")
  TEST_EQUAL(cgp_list[1].params.getValue("recalculate_peaks"), "false")
  TEST_EQUAL(cgp_list[1].params.getValue("use_precursors"), "true")
  TEST_REAL_SIMILAR(cgp_list[1].params.getValue("recalculate_peaks_max_z"), 3.0)
  TEST_REAL_SIMILAR(cgp_list[1].params.getValue("minimal_quality"), -10003.0)
  TEST_REAL_SIMILAR(cgp_list[1].params.getValue("resample_boundary"), 0.03)
  TEST_EQUAL(cgp_list[1].params.getValue("compute_peak_quality"), "true")
  TEST_EQUAL(cgp_list[1].params.getValue("compute_peak_shape_metrics"), "false")

  TEST_EQUAL(cgp_list[3].component_group_name, "ser-L")
  TEST_EQUAL(cgp_list[3].params.getValue("stop_after_feature"), 16)
  TEST_REAL_SIMILAR(cgp_list[3].params.getValue("stop_after_intensity_ratio"), 0.0016)
  TEST_REAL_SIMILAR(cgp_list[3].params.getValue("min_peak_width"), -16.0)
  TEST_EQUAL(cgp_list[3].params.getValue("peak_integration"), "smoothed8")
  TEST_EQUAL(cgp_list[3].params.getValue("background_subtraction"), "none8")
  TEST_EQUAL(cgp_list[3].params.getValue("recalculate_peaks"), "true")
  TEST_EQUAL(cgp_list[3].params.getValue("use_precursors"), "false")
  TEST_REAL_SIMILAR(cgp_list[3].params.getValue("recalculate_peaks_max_z"), 8.0)
  TEST_REAL_SIMILAR(cgp_list[3].params.getValue("minimal_quality"), -10008.0)
  TEST_REAL_SIMILAR(cgp_list[3].params.getValue("resample_boundary"), 0.08)
  TEST_EQUAL(cgp_list[3].params.getValue("compute_peak_quality"), "false")
  TEST_EQUAL(cgp_list[3].params.getValue("compute_peak_shape_metrics"), "true")

  TEST_EQUAL(cgp_list[4].component_group_name, "group2")

  TEST_REAL_SIMILAR(cgp_list[4].params.getValue("stop_after_intensity_ratio"), 0.0026)
  TEST_REAL_SIMILAR(cgp_list[4].params.getValue("min_peak_width"), -26.0)
  TEST_EQUAL(cgp_list[4].params.getValue("peak_integration"), "smoothed13")
  TEST_EQUAL(cgp_list[4].params.getValue("background_subtraction"), "none13")
  TEST_EQUAL(cgp_list[4].params.getValue("recalculate_peaks"), "false")
  TEST_EQUAL(cgp_list[4].params.getValue("use_precursors"), "true")
  TEST_REAL_SIMILAR(cgp_list[4].params.getValue("recalculate_peaks_max_z"), 13.0)
  TEST_REAL_SIMILAR(cgp_list[4].params.getValue("minimal_quality"), -10013.0)
  TEST_REAL_SIMILAR(cgp_list[4].params.getValue("resample_boundary"), 0.13)
  TEST_EQUAL(cgp_list[4].params.getValue("compute_peak_quality"), "true")
  TEST_EQUAL(cgp_list[4].params.getValue("compute_peak_shape_metrics"), "false")

  TEST_EQUAL(cgp_list[4].params.exists("stop_after_feature"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
