// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeaturePicker.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <map>

namespace OpenMS
{
  /**
    @brief Loads MRM component and component-group parameter sets from a comma-separated file.

    Produces parallel @c ComponentParams and @c ComponentGroupParams
    vectors (see @ref MRMFeaturePicker) for use by @ref MRMFeaturePicker.

    @par File format
    Two required columns: @c "component_name" and @c "component_group_name".
    Any additional column is interpreted by its header name:
      - @c "TransitionGroupPicker:PeakPickerChromatogram:&lt;param&gt;" feeds the
        per-component parameter set (@c ComponentParams::params).
      - @c "TransitionGroupPicker:&lt;param&gt;" (without the
        @c PeakPickerChromatogram prefix) feeds the per-component-group
        parameter set (@c ComponentGroupParams::params).
      - Any other column is ignored.
    Cells whose value is empty are not written into the parameter set.
    Lines whose @c component_name or @c component_group_name cell is
    empty are skipped entirely.

    A reduced example (fewer columns are shown here):
    > component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length
    > arg-L.arg-L_1.Heavy,arg-L,2,15
    > arg-L.arg-L_1.Light,arg-L,2,17
    > orn.orn_1.Heavy,orn,3,21
    > orn.orn_1.Light,orn,3,13

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MRMFeaturePickerFile :
    public CsvFile
  {
public:
    /// Default constructor.
    MRMFeaturePickerFile() = default;
    /// Destructor.
    ~MRMFeaturePickerFile() override = default;

    /**
      @brief Read @p filename and populate @p cp_list / @p cgp_list.

      Both output vectors are cleared at the start of the call. The file
      is parsed as comma-separated; the first non-empty line is treated
      as the header. For every subsequent line with non-empty
      @c component_name and @c component_group_name cells one
      @c ComponentParams entry is appended to @p cp_list. A
      @c ComponentGroupParams entry is appended to @p cgp_list only the
      first time a given @c component_group_name is encountered;
      subsequent rows with the same group name extend @p cp_list but do
      not produce a new @p cgp_list entry.

      @note If the file contains zero or one line, no header validation
            is performed and the output vectors stay empty.

      @param[in]  filename Path to the @c .csv input file.
      @param[out] cp_list  Receives the per-row component parameters.
      @param[out] cgp_list Receives one entry per unique @c component_group_name.

      @throws Exception::MissingInformation When the required
                                            @c component_name and/or
                                            @c component_group_name
                                            columns are absent.
      @throws Exception::FileNotFound       When @p filename cannot be
                                            opened.
    */
    void load(
      const String& filename,
      std::vector<MRMFeaturePicker::ComponentParams>& cp_list,
      std::vector<MRMFeaturePicker::ComponentGroupParams>& cgp_list
    );

protected:
    /**
      @brief Extract the @c ComponentParams and @c ComponentGroupParams entries from a single parsed CSV line.

      Reads @c component_name and @c component_group_name from @p line via
      @p headers; if either is empty, the call returns @c false and
      @p cp / @p cgp are left in an unspecified state. Otherwise every
      column whose header matches one of the @c TransitionGroupPicker:
      prefixes feeds the corresponding parameter set as described in the
      class brief.

      @param[in]  line    Tokens of the input line, indexed by @p headers.
      @param[in]  headers Mapping from header name to column index in @p line.
      @param[out] cp      Receives the per-component parameters.
      @param[out] cgp     Receives the per-component-group parameters (@c component_group_name plus any matching values).

      @return @c true on a fully-populated pair; @c false when either
              required name cell was empty.
    */
    bool extractParamsFromLine_(
      const StringList& line,
      const std::map<String, Size>& headers,
      MRMFeaturePicker::ComponentParams& cp,
      MRMFeaturePicker::ComponentGroupParams& cgp
    ) const;

    /**
      @brief Insert @p value into @p params under @p key, converting the string to the type expected for that parameter name.

      Empty values are dropped without modifying @p params. The target
      type is selected by looking @p key up in a hard-coded table:
        - @c "gauss_width", @c "peak_width", @c "signal_to_noise",
          @c "sn_win_len", @c "stop_after_intensity_ratio",
          @c "min_peak_width", @c "recalculate_peaks_max_z",
          @c "minimal_quality", @c "resample_boundary" -> @c double.
        - @c "use_gauss", @c "write_sn_log_messages",
          @c "remove_overlapping_peaks", @c "recalculate_peaks",
          @c "use_precursors", @c "compute_peak_quality",
          @c "compute_peak_shape_metrics" -> @c bool (case-insensitive
          @c "true" maps to @c "true", everything else to @c "false").
        - @c "sgolay_frame_length", @c "sgolay_polynomial_order",
          @c "sn_bin_count" -> @c UInt (via @c double rounding).
        - @c "stop_after_feature" -> @c Int.
        - Anything else -> stored as a @c String.

      @param[in]     key    Parameter name (the column header with any @c TransitionGroupPicker: prefix already stripped).
      @param[in]     value  Value to convert.
      @param[in,out] params Parameter container to update.
    */
    void setCastValue_(const String& key, const String& value, Param& params) const;
  };
}

