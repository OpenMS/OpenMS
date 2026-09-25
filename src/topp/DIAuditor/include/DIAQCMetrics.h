// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, David L. Tabb $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <array>
#include <iosfwd>
#include <string>
#include <vector>

namespace OpenMS
{
  class ExperimentalSettings;
  class MSSpectrum;

  /**
    @brief Quality metrics for data-independent acquisition (DIA) runs, per run and per isolation window.

    Re-implements the metrics of DIAuditor by David L. Tabb (https://github.com/dtabb73/DIAuditor) on top of the
    OpenMS mzML reader. Only a small record per spectrum is kept, so a run can be streamed:

    @code
    DIAQCMetrics metrics;
    MSDataTransformingConsumer consumer;
    consumer.setExperimentalSettingsFunc([&](const ExperimentalSettings& s) { metrics.setExperimentalSettings(s); });
    consumer.setSpectraProcessingFunc([&](MSSpectrum& s) { metrics.addSpectrum(s); });
    MzMLFile().transform(file, &consumer, true);
    DIAQCMetrics::RunMetrics run = metrics.compute();
    @endcode

    MS2 spectra are grouped into isolation windows. Two spectra belong to the same window if their first precursor has
    the same isolation window (target m/z, lower and upper offset) and, depending on Options::ion_mobility, the same
    FAIMS compensation voltage and ion mobility range. Values are compared with an absolute tolerance of 1e-6 (as in
    OpenSWATH), since converters may write them with different last digits. Windows are reported in the order in which
    they are first acquired. Spectra of MS level 3 and higher are counted, but not assigned to windows. Spectra without
    a retention time (no scan start time) are counted, but not used for any other metric.

    Conventions (retention times are kept in seconds; the table writers convert them to minutes like DIAuditor):
    - Medians of times are the usual median (mean of the two middle values for an even count).
    - Peak count quartiles, including their median, are observed values taken as in DIAuditor: the values at
      positions n/4, n/2 and n/4 + n/2 (integer division) of the n sorted counts.
    - A "TIC quantile RT" is the retention time of the first spectrum at which the cumulative TIC reaches the quantile.
    - Cycle times are medians of the time between consecutive spectra of the same kind (MS1, or one isolation window).
    - Undefined values (e.g. a cycle time from a single spectrum, TIC quantiles when the TIC is zero) are NaN.
  */
  class DIAQCMetrics
  {
  public:
    /// Where the total ion current of a spectrum comes from
    enum class TICSource
    {
      AUTO,     ///< the 'total ion current' value of the file (MS:1000285) if present, otherwise the sum of intensities
      FILE,     ///< only the 'total ion current' value of the file; spectra without it count with a TIC of 0 (like DIAuditor)
      COMPUTED  ///< always the sum of intensities of the spectrum
    };

    /// How peaks of a spectrum are counted
    enum class PeakCountMode
    {
      ALL,      ///< all data points (like DIAuditor's defaultArrayLength), including zero-intensity points of profile data
      NONZERO   ///< only data points with an intensity above zero
    };

    /// Which ion mobility separates isolation windows with the same m/z range
    enum class IonMobilityKey
    {
      AUTO,     ///< FAIMS compensation voltage and the ion mobility range of the window ('ion mobility lower/upper limit' of
                ///< e.g. diaPASEF frames). An ion mobility of a single spectrum (e.g. one TIMS scan) is a position, not a window.
      FAIMS,    ///< only the FAIMS compensation voltage (like DIAuditor)
      NONE      ///< none: windows are defined by their m/z range alone
    };

    /// Options
    struct Options
    {
      TICSource tic_source = TICSource::AUTO;
      PeakCountMode peak_count = PeakCountMode::ALL;
      IonMobilityKey ion_mobility = IonMobilityKey::AUTO;
    };

    /// Minimum, first quartile, median, third quartile and maximum of peak counts
    using PeakCountSummary = std::array<double, 5>;

    /// Retention times (seconds) at which 25%, 50% and 75% of the TIC have been accumulated
    using TICQuantileRTs = std::array<double, 3>;

    /// Metrics of one isolation window
    struct WindowMetrics
    {
      bool has_isolation_window = false;  ///< false for MS2 spectra without isolation window offsets
      double target_mz = 0.0;             ///< isolation window target (or the precursor m/z if there is no isolation window)
      double lower_mz = 0.0;              ///< target m/z minus lower offset
      double upper_mz = 0.0;              ///< target m/z plus upper offset
      double width_mz = 0.0;              ///< upper_mz - lower_mz
      double faims_cv;                    ///< FAIMS compensation voltage; NaN if none
      double ion_mobility_lower;          ///< lower limit of the ion mobility range of the window; NaN if none
      double ion_mobility_upper;          ///< upper limit of the ion mobility range of the window; NaN if none
      double mass_resolving_power;        ///< median over the spectra of the window; NaN if not annotated
      Size spectrum_count = 0;
      double rt_min;                      ///< seconds
      double rt_max;                      ///< seconds
      double cycle_time_median;           ///< seconds
      TICQuantileRTs tic_quantile_rt;
      double total_tic = 0.0;
      PeakCountSummary peak_count;

      WindowMetrics();
    };

    /// Metrics of one run, including its isolation windows
    struct RunMetrics
    {
      std::string source_file;            ///< file name without directory and extension
      std::string input_path;             ///< path of the input file (for mzQC)
      std::string file_sha1;              ///< SHA-1 of the input file (for mzQC); may be empty
      std::string instrument;             ///< instrument model name
      std::string serial_number;          ///< instrument serial number
      std::string start_time_stamp;       ///< start of the acquisition as written in the file

      double rt_min;                      ///< seconds, over all spectra
      double rt_max;                      ///< seconds, over all spectra
      Size ms1_count = 0;
      Size msn_count = 0;                 ///< spectra of MS level 2 and higher
      Size ms2_count = 0;                 ///< spectra of MS level 2
      Size spectra_without_rt = 0;        ///< spectra without retention time (not used for other metrics)
      Size ms2_multiple_precursors = 0;   ///< MS2 spectra with more than one precursor (e.g. multiplexed DIA); only the first is used
      Size ms2_scan_ion_mobility = 0;     ///< MS2 spectra that look like single ion mobility scans (no range or IM array)

      double ms1_mass_resolving_power;    ///< median over MS1 spectra
      TICQuantileRTs ms1_tic_quantile_rt;
      double ms1_total_tic = 0.0;
      double ms1_cycle_time_median;
      PeakCountSummary ms1_peak_count;

      TICQuantileRTs ms2_tic_quantile_rt;
      double ms2_total_tic = 0.0;
      PeakCountSummary ms2_peak_count;

      Size window_count = 0;
      Size windows_measured_once = 0;     ///< windows with a single MS2 spectrum
      double window_spectra_min;          ///< fewest spectra of any window
      double window_spectra_max;          ///< most spectra of any window
      double window_mz_min;               ///< lowest lower m/z of any window
      double window_mz_max;               ///< highest upper m/z of any window
      double window_width_min;
      double window_width_max;
      double window_cycle_time_mean;      ///< mean of the windows' median cycle times (DIAuditor's AverageMedianCycleTime)
      double window_cycle_time_median;    ///< median of all times between consecutive spectra of the same window
      double window_half_tic_rt_min;      ///< earliest of the windows' 50% TIC retention times
      double window_half_tic_rt_max;      ///< latest of the windows' 50% TIC retention times
      double window_total_tic_min;
      double window_total_tic_max;
      double window_peak_count_median_min;
      double window_peak_count_median_max;

      std::vector<WindowMetrics> windows;

      RunMetrics();
    };

    /// Default options
    DIAQCMetrics();

    explicit DIAQCMetrics(const Options& options);

    /// Take instrument, serial number and start time from the run's settings
    void setExperimentalSettings(const ExperimentalSettings& settings);

    /// Record one spectrum; spectra can come in any order
    void addSpectrum(const MSSpectrum& spectrum);

    /// Number of spectra recorded so far
    Size size() const;

    /// Compute the metrics of the recorded run
    RunMetrics compute() const;

    /// Forget all recorded spectra and settings
    void clear();

    /// Write one row per run (DIAuditor's "byRun" table) as tab-separated values
    static void writeRunTable(const std::vector<RunMetrics>& runs, std::ostream& os);

    /// Write one row per isolation window (DIAuditor's "byIsolationWindow" table) as tab-separated values
    static void writeWindowTable(const std::vector<RunMetrics>& runs, std::ostream& os);

    /**
      @brief Write an mzQC 1.0 file with one runQuality per run

      Only metrics defined in the PSI-MS vocabulary are written, with the names and units the vocabulary defines. Terms
      missing from the installed vocabulary are skipped with a warning. Per-window values are not written, as the
      vocabulary has no term for them; see writeWindowTable().

      @param[in] runs The runs
      @param[in] os Output stream
      @param[in] software_version Version string of the software written as analysisSoftware
      @param[in] creation_date Creation date as RFC 3339 date-time, with time zone (e.g. "2026-01-02T03:04:05Z")
    */
    static void writeMzQC(const std::vector<RunMetrics>& runs, std::ostream& os, const std::string& software_version, const std::string& creation_date);

  private:
    /// The per-spectrum data the metrics need
    struct SpectrumRecord
    {
      double rt = 0.0;
      UInt ms_level = 0;
      double tic = 0.0;
      Size peak_count = 0;
      double mass_resolving_power = 0.0;  // NaN if not annotated
      bool has_isolation_window = false;
      double target_mz = 0.0;
      double lower_offset = 0.0;
      double upper_offset = 0.0;
      double faims_cv = 0.0;              // NaN if none
      double ion_mobility_lower = 0.0;    // NaN if none
      double ion_mobility_upper = 0.0;    // NaN if none
      bool scan_ion_mobility = false;     // the spectrum looks like a single ion mobility scan
      Size precursor_count = 0;
    };

    /// Statistics of a group of records (all MS1, all MS2, or one isolation window); defined in the source file
    struct GroupStatistics;

    Options options_;
    std::vector<SpectrumRecord> records_;
    std::string instrument_;
    std::string serial_number_;
    std::string start_time_stamp_;
  };
} // namespace OpenMS
