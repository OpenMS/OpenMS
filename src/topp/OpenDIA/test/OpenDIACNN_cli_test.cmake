# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Justin Sing $
# $Authors: Justin Sing $
# --------------------------------------------------------------------------

file(MAKE_DIRECTORY "${WORK}/input.d")
file(COPY_FILE "${LIBRARY}" "${WORK}/library.pqp")
file(WRITE "${WORK}/input.d/analysis.tdf" "raw input sentinel\n")
file(WRITE "${WORK}/windows.tsv" "sample_id\tprecursor_id\tgroup_id\tlabel\trt_apex\n")
foreach(model IN ITEMS "${WORK}/library.pqp" "${WORK}/input.d/analysis.tdf")
  file(SHA256 "${model}" before)
  execute_process(COMMAND "${BINARY}" -mode train -in "${WORK}/input.d"
    -lib "${WORK}/library.pqp" -windows "${WORK}/windows.tsv"
    -model "${model}" -out "${WORK}/scores.tsv"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "must not overwrite")
    message(FATAL_ERROR "Native CLI did not reject input overwrite: ${status}\n${output}\n${error}")
  endif()
  file(SHA256 "${model}" after)
  if(NOT before STREQUAL after OR EXISTS "${WORK}/scores.tsv")
    message(FATAL_ERROR "Native CLI changed input or opened score output before validation")
  endif()
endforeach()

# Crossfit errors must be detected before touching raw data or score outputs.
file(WRITE "${WORK}/all.tsv" "sample_id\tprecursor_id\tgroup_id\tlabel\trt_apex\none\t1\tfamily\t-1\t10\n")
file(WRITE "${WORK}/seeds.tsv" "sample_id\tprecursor_id\tgroup_id\tlabel\trt_apex\none\t1\tfamily\t1\t11\n")
execute_process(COMMAND "${BINARY}" -mode crossfit -in "${WORK}/input.d"
  -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -seeds "${WORK}/seeds.tsv"
  -model "${WORK}/crossfit" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "RT coordinates must match exactly" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "Crossfit did not reject mismatched seed geometry before extraction: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -in "${WORK}/input.d"
  -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv"
  -model "${WORK}/crossfit" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "Crossfit requires -seeds" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "Crossfit did not require an explicit seed manifest: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -in "${WORK}/input.d"
  -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -seeds "${WORK}/seeds.tsv"
  -model "${WORK}/input.d/checkpoint" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "must not overwrite" OR EXISTS "${WORK}/input.d/checkpoint.fold0.pt")
  message(FATAL_ERROR "Crossfit did not protect the raw directory from derived checkpoint paths")
endif()
file(WRITE "${WORK}/existing.fold0.pt" "checkpoint sentinel\n")
file(SHA256 "${WORK}/existing.fold0.pt" before)
execute_process(COMMAND "${BINARY}" -mode crossfit -in "${WORK}/input.d"
  -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -seeds "${WORK}/seeds.tsv"
  -model "${WORK}/existing" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
file(SHA256 "${WORK}/existing.fold0.pt" after)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "refuses to overwrite" OR NOT before STREQUAL after)
  message(FATAL_ERROR "Crossfit overwrote or failed to reject an existing fold checkpoint")
endif()

# Paired competition is crossfit-only and balances the classes itself (ratio must stay 1).
foreach(arguments IN ITEMS "-mode;train;-model;${WORK}/paired.pt" "-mode;crossfit;-seeds;${WORK}/seeds.tsv;-model;${WORK}/paired;-target_decoy_ratio;1.5")
  execute_process(COMMAND "${BINARY}" ${arguments} -paired_tdc -in "${WORK}/input.d"
    -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "paired_tdc requires crossfit" OR EXISTS "${WORK}/scores.tsv")
    message(FATAL_ERROR "paired_tdc misuse was not rejected before extraction: ${arguments}\n${output}\n${error}")
  endif()
endforeach()

# Pooled competition only refines paired competition.
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/pooled" -pooled_tdc
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "pooled_tdc requires paired_tdc" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "pooled_tdc without paired_tdc was not rejected before extraction\n${output}\n${error}")
endif()

# The experimental RT window is crossfit-only and is rejected before extraction.
execute_process(COMMAND "${BINARY}" -mode train -rt_window_calibration -in "${WORK}/input.d"
  -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -model "${WORK}/rt.pt" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "rt_window_calibration requires crossfit" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "rt_window_calibration outside crossfit was not rejected: ${status}\n${output}\n${error}")
endif()

# The identification-free calibration file is crossfit-only, excludes the CNN-anchor window, and an invalid
# file fails before extraction.
file(WRITE "${WORK}/calibration.json" "{\"version\": 1, \"rt\": {\"knots_library_rt\": [0], \"knots_run_rt_s\": [10], \"sigma_s\": 1, \"anchors_used\": 20, \"tolerance_half_width_s\": 4}}\n")
foreach(arguments IN ITEMS "-mode;train;-model;${WORK}/cal.pt" "-mode;crossfit;-seeds;${WORK}/seeds.tsv;-model;${WORK}/cal;-rt_window_calibration")
  execute_process(COMMAND "${BINARY}" ${arguments} -rt_calibration "${WORK}/calibration.json" -in "${WORK}/input.d"
    -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "rt_calibration requires crossfit" OR EXISTS "${WORK}/scores.tsv")
    message(FATAL_ERROR "rt_calibration misuse was not rejected before extraction: ${arguments}\n${output}\n${error}")
  endif()
endforeach()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/cal" -rt_calibration "${WORK}/calibration.json"
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "invalid RT map" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "An invalid calibration file was not rejected before extraction: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/cal" -apply_calibration_coordinates
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "apply_calibration_coordinates requires -rt_calibration" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "apply_calibration_coordinates without a calibration was not rejected: ${status}\n${output}\n${error}")
endif()

# Score iterations are crossfit-only and are rejected before extraction.
execute_process(COMMAND "${BINARY}" -mode train -ss_iterations 1 -model "${WORK}/ss.pt"
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "ss_iterations requires crossfit" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "ss_iterations outside crossfit was not rejected: ${status}\n${output}\n${error}")
endif()

# The second-stage score needs pooled paired competition, and its feature table needs -stage2; both fail before extraction.
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/stage2" -paired_tdc -stage2
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "stage2 requires crossfit mode with -paired_tdc and -pooled_tdc" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "stage2 without pooled_tdc was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/stage2f" -paired_tdc -pooled_tdc
  -stage2_features_out "${WORK}/features.tsv"
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "stage2_features_out requires -stage2" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "stage2_features_out without stage2 was not rejected: ${status}\n${output}\n${error}")
endif()
# Search peak-group features need -stage2 and a matching feature set; both fail before extraction.
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/stage2s" -paired_tdc -pooled_tdc -stage2
  -stage2_feature_set search
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "stage2_search_features is required by" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "stage2_feature_set search without a search table was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/stage2t" -paired_tdc -pooled_tdc
  -stage2_search_features "${WORK}/all.tsv"
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "require -stage2" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "stage2_search_features without stage2 was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/stage2u" -paired_tdc -pooled_tdc -stage2
  -stage2_feature_set search -stage2_search_features "${WORK}/all.tsv"
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "exactly one column candidate_row" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "malformed search table was not rejected before extraction: ${status}\n${output}\n${error}")
endif()
# Multiple-instance training needs -cnn_arch topaz and an identification-free calibration; both fail before extraction.
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/mil" -cnn_training mil
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "cnn_training mil requires" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "cnn_training mil without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/mil2" -cnn_arch topaz -cnn_training mil
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "cnn_training mil requires" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "cnn_training mil without rt_calibration was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/arch" -cnn_arch resnet
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "unknown cnn_arch was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/width" -topaz_width 8
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "topaz_width, weight_decay and topaz_faithful require -cnn_arch topaz" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_width without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/faithful" -topaz_faithful
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "topaz_faithful require -cnn_arch topaz" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_faithful without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
# Single-scorer TOPAZ heuristic features and augmentation fail before extraction when misused.
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/tr1" -trace_refinement_min_r 0.75
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "trace_refinement_min_r requires -stage2" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "trace_refinement_min_r without stage2/topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/tf1" -topaz_features byion
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "topaz_features requires -cnn_arch topaz" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_features without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/tf2" -cnn_arch topaz -topaz_features byion
  -paired_tdc -pooled_tdc -stage2
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "excludes -stage2" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_features with stage2 was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/tf3" -cnn_arch topaz -topaz_features search
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "stage2_search_features is required by" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_features search without a search table was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/tf4" -cnn_arch topaz -topaz_features nonsense
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "Unknown or repeated" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "unknown topaz_features group was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/aug1" -augment_probability 0.5
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "augment_probability requires -cnn_arch topaz" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "augment_probability without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/dn1" -topaz_decoy_negatives 100
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "topaz_decoy_negatives requires -cnn_arch topaz" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_decoy_negatives without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/lh1" -topaz_linear_head
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "topaz_linear_head requires -cnn_arch topaz" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_linear_head without topaz was not rejected: ${status}\n${output}\n${error}")
endif()
execute_process(COMMAND "${BINARY}" -mode crossfit -seeds "${WORK}/seeds.tsv" -model "${WORK}/ns1" -cnn_arch topaz -topaz_negative_source target
  -in "${WORK}/input.d" -lib "${WORK}/library.pqp" -windows "${WORK}/all.tsv" -out "${WORK}/scores.tsv"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status STREQUAL "0" OR NOT "${output}${error}" MATCHES "topaz_negative_source requires -topaz_decoy_negatives" OR EXISTS "${WORK}/scores.tsv")
  message(FATAL_ERROR "topaz_negative_source without decoy negatives was not rejected: ${status}\n${output}\n${error}")
endif()
