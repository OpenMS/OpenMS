// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, David L. Tabb $
// --------------------------------------------------------------------------

#include "DIAQCMetrics.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <utility>

namespace OpenMS
{
  namespace
  {
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
    constexpr std::array<double, 3> TIC_QUANTILES = {0.25, 0.5, 0.75};
    /// absolute tolerance for values that define an isolation window (as in FullSwathFileConsumer)
    constexpr double WINDOW_TOLERANCE = 1e-6;

    /// A numeric meta value as double; NaN if missing or not numeric
    double metaValueAsDouble(const MetaInfoInterface& meta, const std::string& name)
    {
      if (!meta.metaValueExists(name)) return NaN;
      const DataValue& value = meta.getMetaValue(name);
      switch (value.valueType())
      {
        case DataValue::DOUBLE_VALUE:
        case DataValue::INT_VALUE:
          return static_cast<double>(value);
        case DataValue::STRING_VALUE:
          try
          {
            return StringUtils::toDouble(value.toString());
          }
          catch (const Exception::ConversionError&)
          {
            return NaN;
          }
        default:
          return NaN;
      }
    }

    /// Minimum, quartiles and maximum of sorted counts, taken as in DIAuditor: the values at positions n/4, n/2 and
    /// n/4 + n/2 (integer division)
    DIAQCMetrics::PeakCountSummary peakCountSummary(const std::vector<Size>& sorted)
    {
      if (sorted.empty()) return {NaN, NaN, NaN, NaN, NaN};
      const Size n = sorted.size();
      return {static_cast<double>(sorted.front()), static_cast<double>(sorted[n / 4]), static_cast<double>(sorted[n / 2]),
              static_cast<double>(sorted[n / 4 + n / 2]), static_cast<double>(sorted.back())};
    }

    /// Sum of intensities, accumulated in double precision
    double intensitySum(const MSSpectrum& spectrum)
    {
      double sum = 0.0;
      for (const Peak1D& p : spectrum) sum += p.getIntensity();
      return sum;
    }

    /// Two values of a window definition are the same (both undefined, or equal within the tolerance)
    bool sameValue(double a, double b)
    {
      return (std::isnan(a) && std::isnan(b)) || std::fabs(a - b) <= WINDOW_TOLERANCE;
    }

    double medianOrNaN(std::vector<double> values)
    {
      values.erase(std::remove_if(values.begin(), values.end(), [](double v) { return !std::isfinite(v); }), values.end());
      if (values.empty()) return NaN;
      return Math::median(values.begin(), values.end());
    }

    /// Minimum and maximum of the finite values; NaN if there are none
    std::pair<double, double> finiteRange(const std::vector<double>& values)
    {
      double lo = NaN, hi = NaN;
      for (double v : values)
      {
        if (!std::isfinite(v)) continue;
        if (!(v >= lo)) lo = v; // also true while lo is NaN
        if (!(v <= hi)) hi = v;
      }
      return {lo, hi};
    }

    /// Up to 12 significant digits, independent of the locale: plain notation for the usual magnitudes (120000, not
    /// "1.2e05") and no noise from double arithmetic (0.06, not 0.0600000000000094)
    std::string numberToString(double value)
    {
      if (!std::isfinite(value)) return "NA";
      if (value == 0.0) value = 0.0; // no "-0"
      std::array<char, 32> buffer;
      const auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), value, std::chars_format::general, 12);
      return std::string(buffer.data(), result.ptr);
    }

    std::string countToString(double value)
    {
      if (!std::isfinite(value)) return "NA";
      return StringUtils::toStr(static_cast<long long>(std::llround(value)));
    }

    double toMinutes(double seconds)
    {
      return seconds / 60.0;
    }

    /// Text for a table cell: tabs and line breaks would break the table
    std::string textCell(const std::string& text)
    {
      if (text.empty()) return "NA";
      std::string cell(text);
      std::replace_if(cell.begin(), cell.end(), [](char c) { return c == '\t' || c == '\n' || c == '\r'; }, ' ');
      return cell;
    }

    void writeRow(std::ostream& os, const std::vector<std::string>& cells)
    {
      for (Size i = 0; i < cells.size(); ++i)
      {
        if (i > 0) os << '\t';
        os << cells[i];
      }
      os << '\n';
    }

    /// data-version of an OBO file (empty if not found); reads the header only
    std::string oboDataVersion(const std::string& obo_file)
    {
      std::ifstream is(obo_file);
      std::string line;
      const std::string key = "data-version:";
      while (std::getline(is, line))
      {
        if (line.starts_with(key)) return StringUtils::trimmed(line.substr(key.size()));
        if (line.starts_with("[Term]")) break; // end of the header
      }
      return "";
    }

    /// PSI-MS term (accession, name) of an input file format, for the mzQC inputFile
    std::pair<std::string, std::string> fileFormatTerm(FileTypes::Type type)
    {
      switch (type)
      {
        case FileTypes::RAW: return {"MS:1000563", "Thermo RAW format"};
        case FileTypes::BRUKER_TDF: return {"MS:1002817", "Bruker TDF format"};
        default: return {"MS:1000584", "mzML format"};
      }
    }

    /// file URI of a local path: every byte of its UTF-8 form outside the unreserved characters of RFC 3986 and '/' is
    /// percent-encoded (a backslash is a separator only on Windows, where generic_u8string() turns it into '/')
    std::string fileURI(const std::string& path)
    {
      const std::u8string utf8 = std::filesystem::absolute(to_path(path)).lexically_normal().generic_u8string();
      std::string absolute(reinterpret_cast<const char*>(utf8.data()), utf8.size());
      if (!absolute.starts_with("/")) absolute.insert(0, "/"); // Windows drive letter
      const bool drive = absolute.size() > 2 && std::isalpha(static_cast<unsigned char>(absolute[1])) && absolute[2] == ':';
      std::string uri = "file://";
      const char* hex = "0123456789ABCDEF";
      for (Size i = 0; i < absolute.size(); ++i)
      {
        const unsigned char c = static_cast<unsigned char>(absolute[i]);
        if (std::isalnum(c) || c == '-' || c == '.' || c == '_' || c == '~' || c == '/' || (drive && i == 2))
        {
          uri += static_cast<char>(c);
        }
        else
        {
          uri += '%';
          uri += hex[c >> 4];
          uri += hex[c & 0x0F];
        }
      }
      return uri;
    }
  } // namespace

  struct DIAQCMetrics::GroupStatistics
  {
    double total_tic = 0.0;
    TICQuantileRTs tic_quantile_rt{NaN, NaN, NaN};
    PeakCountSummary peak_count{NaN, NaN, NaN, NaN, NaN};
    std::vector<double> rt_differences;
    double cycle_time_median = NaN;
    double mass_resolving_power = NaN;
    double rt_min = NaN;
    double rt_max = NaN;

    /// @p records must be sorted by retention time
    explicit GroupStatistics(const std::vector<const SpectrumRecord*>& records)
    {
      if (records.empty()) return;
      rt_min = records.front()->rt;
      rt_max = records.back()->rt;

      for (const SpectrumRecord* r : records) total_tic += r->tic;
      if (total_tic > 0.0)
      {
        double cumulative = 0.0;
        Size next = 0;
        for (const SpectrumRecord* r : records)
        {
          cumulative += r->tic;
          while (next < TIC_QUANTILES.size() && cumulative >= TIC_QUANTILES[next] * total_tic)
          {
            tic_quantile_rt[next++] = r->rt;
          }
        }
      }

      std::vector<Size> counts;
      counts.reserve(records.size());
      for (const SpectrumRecord* r : records) counts.push_back(r->peak_count);
      std::sort(counts.begin(), counts.end());
      peak_count = peakCountSummary(counts);

      rt_differences.reserve(records.size());
      for (Size i = 1; i < records.size(); ++i) rt_differences.push_back(records[i]->rt - records[i - 1]->rt);
      cycle_time_median = medianOrNaN(rt_differences);

      std::vector<double> resolving_powers;
      for (const SpectrumRecord* r : records) resolving_powers.push_back(r->mass_resolving_power);
      mass_resolving_power = medianOrNaN(std::move(resolving_powers));
    }
  };

  DIAQCMetrics::WindowMetrics::WindowMetrics() :
    target_mz(NaN),
    lower_mz(NaN),
    upper_mz(NaN),
    width_mz(NaN),
    faims_cv(NaN),
    ion_mobility_lower(NaN),
    ion_mobility_upper(NaN),
    mass_resolving_power(NaN),
    rt_min(NaN),
    rt_max(NaN),
    cycle_time_median(NaN),
    tic_quantile_rt{NaN, NaN, NaN},
    peak_count{NaN, NaN, NaN, NaN, NaN}
  {
  }

  DIAQCMetrics::RunMetrics::RunMetrics() :
    rt_min(NaN),
    rt_max(NaN),
    precursor_mz_min(NaN),
    precursor_mz_max(NaN),
    ms1_mass_resolving_power(NaN),
    ms1_tic_quantile_rt{NaN, NaN, NaN},
    ms1_cycle_time_median(NaN),
    ms1_peak_count{NaN, NaN, NaN, NaN, NaN},
    ms2_tic_quantile_rt{NaN, NaN, NaN},
    ms2_peak_count{NaN, NaN, NaN, NaN, NaN},
    window_spectra_min(NaN),
    window_spectra_max(NaN),
    window_mz_min(NaN),
    window_mz_max(NaN),
    window_width_min(NaN),
    window_width_max(NaN),
    window_cycle_time_mean(NaN),
    window_cycle_time_median(NaN),
    window_half_tic_rt_min(NaN),
    window_half_tic_rt_max(NaN),
    window_total_tic_min(NaN),
    window_total_tic_max(NaN),
    window_peak_count_median_min(NaN),
    window_peak_count_median_max(NaN)
  {
  }

  DIAQCMetrics::DIAQCMetrics() = default;

  DIAQCMetrics::DIAQCMetrics(const Options& options) :
    options_(options)
  {
  }

  void DIAQCMetrics::setExperimentalSettings(const ExperimentalSettings& settings)
  {
    instrument_ = settings.getInstrument().getName();
    serial_number_.clear();
    if (settings.getInstrument().metaValueExists("instrument serial number"))
    {
      serial_number_ = settings.getInstrument().getMetaValue("instrument serial number").toString();
    }
    // The mzML reader keeps the attribute verbatim if it carries more than DateTime can hold (e.g. a time zone).
    start_time_stamp_.clear();
    if (settings.metaValueExists("mzml_start_time_stamp"))
    {
      start_time_stamp_ = settings.getMetaValue("mzml_start_time_stamp").toString();
    }
    else if (!settings.getDateTime().isNull() && settings.getDateTime().isValid())
    {
      start_time_stamp_ = settings.getDateTime().toString("yyyy-MM-ddThh:mm:ss");
    }
  }

  void DIAQCMetrics::addSpectrum(const MSSpectrum& spectrum)
  {
    SpectrumRecord record;
    record.rt = spectrum.getRT();
    record.ms_level = spectrum.getMSLevel();

    const double file_tic = metaValueAsDouble(spectrum, "total ion current");
    switch (options_.tic_source)
    {
      case TICSource::AUTO:
        record.tic = std::isfinite(file_tic) ? file_tic : intensitySum(spectrum);
        break;
      case TICSource::FILE:
        record.tic = std::isfinite(file_tic) ? file_tic : 0.0;
        break;
      case TICSource::COMPUTED:
        record.tic = intensitySum(spectrum);
        break;
    }

    if (options_.peak_count == PeakCountMode::ALL)
    {
      record.peak_count = spectrum.size();
    }
    else
    {
      record.peak_count = std::count_if(spectrum.begin(), spectrum.end(), [](const Peak1D& p) { return p.getIntensity() > 0; });
    }

    record.mass_resolving_power = metaValueAsDouble(spectrum, "mass resolving power");

    record.precursor_count = spectrum.getPrecursors().size();
    record.precursor_mz = NaN;
    if (record.ms_level >= 2 && record.precursor_count > 0)
    {
      const Precursor& precursor = spectrum.getPrecursors().front();
      // the recorded precursor m/z: the selected ion m/z, or the isolation window target if there is no selected ion
      if (precursor.getMZ() > 0.0) record.precursor_mz = precursor.getMZ();
      // if the file has a selected ion m/z, the mzML reader moves the isolation window target into a meta value
      // (0 if the precursor has no isolation window)
      const double target = metaValueAsDouble(precursor, "isolation window target m/z");
      record.target_mz = (target > 0.0) ? target : precursor.getMZ();
      record.lower_offset = precursor.getIsolationWindowLowerOffset();
      record.upper_offset = precursor.getIsolationWindowUpperOffset();
      record.has_isolation_window = record.lower_offset > 0.0 || record.upper_offset > 0.0;
    }

    // FAIMS: the compensation voltage is a setting of the window
    const bool faims = spectrum.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
    record.faims_cv = (faims && options_.ion_mobility != IonMobilityKey::NONE) ? spectrum.getDriftTime() : NaN;
    // other ion mobility: the range of the window (e.g. a diaPASEF frame) is a setting of the window, while the ion
    // mobility of a single spectrum (e.g. one TIMS scan) is a position within it
    const double im_lower = metaValueAsDouble(spectrum, "ion mobility lower limit");
    const double im_upper = metaValueAsDouble(spectrum, "ion mobility upper limit");
    const bool im_range = std::isfinite(im_lower) || std::isfinite(im_upper);
    record.ion_mobility_lower = (options_.ion_mobility == IonMobilityKey::AUTO) ? im_lower : NaN;
    record.ion_mobility_upper = (options_.ion_mobility == IonMobilityKey::AUTO) ? im_upper : NaN;
    // a single scan has an ion mobility of its own, but neither the range nor the ion mobility array of a frame
    record.scan_ion_mobility = !faims && spectrum.getDriftTimeUnit() != DriftTimeUnit::NONE &&
                               spectrum.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET && !im_range &&
                               !spectrum.containsIMData();

    records_.push_back(record);
  }

  Size DIAQCMetrics::size() const
  {
    return records_.size();
  }

  void DIAQCMetrics::clear()
  {
    records_.clear();
    instrument_.clear();
    serial_number_.clear();
    start_time_stamp_.clear();
  }

  DIAQCMetrics::RunMetrics DIAQCMetrics::compute() const
  {
    RunMetrics run;
    run.instrument = instrument_;
    run.serial_number = serial_number_;
    run.start_time_stamp = start_time_stamp_;

    std::vector<const SpectrumRecord*> ms1, ms2;
    for (const SpectrumRecord& r : records_)
    {
      if (r.ms_level == 1) ++run.ms1_count;
      if (r.ms_level >= 2) ++run.msn_count;
      if (r.ms_level == 2)
      {
        ++run.ms2_count;
        if (r.precursor_count > 1) ++run.ms2_multiple_precursors;
        if (r.scan_ion_mobility) ++run.ms2_scan_ion_mobility;
      }
      // a spectrum without scan start time has a negative retention time
      if (!(r.rt >= 0.0))
      {
        ++run.spectra_without_rt;
        continue;
      }
      if (!(r.rt >= run.rt_min)) run.rt_min = r.rt;
      if (!(r.rt <= run.rt_max)) run.rt_max = r.rt;
      if (r.ms_level >= 2 && std::isfinite(r.precursor_mz))
      {
        if (!(r.precursor_mz >= run.precursor_mz_min)) run.precursor_mz_min = r.precursor_mz;
        if (!(r.precursor_mz <= run.precursor_mz_max)) run.precursor_mz_max = r.precursor_mz;
      }
      if (r.ms_level == 1)
      {
        ms1.push_back(&r);
      }
      else if (r.ms_level == 2)
      {
        ms2.push_back(&r);
      }
    }
    auto by_rt = [](const SpectrumRecord* a, const SpectrumRecord* b) { return a->rt < b->rt; };
    std::stable_sort(ms1.begin(), ms1.end(), by_rt);
    std::stable_sort(ms2.begin(), ms2.end(), by_rt);

    const GroupStatistics ms1_stats(ms1);
    run.ms1_mass_resolving_power = ms1_stats.mass_resolving_power;
    run.ms1_tic_quantile_rt = ms1_stats.tic_quantile_rt;
    run.ms1_total_tic = ms1_stats.total_tic;
    run.ms1_cycle_time_median = ms1_stats.cycle_time_median;
    run.ms1_peak_count = ms1_stats.peak_count;

    const GroupStatistics ms2_stats(ms2);
    run.ms2_tic_quantile_rt = ms2_stats.tic_quantile_rt;
    run.ms2_total_tic = ms2_stats.total_tic;
    run.ms2_peak_count = ms2_stats.peak_count;

    // group MS2 spectra into isolation windows, in order of first acquisition; the values that define a window are
    // matched within the tolerance, and candidate windows are looked up by target m/z in bins of 0.001 (much wider).
    // MS2 spectra without an isolation window are not a DIA isolation window and are collected separately.
    auto sameWindow = [](const SpectrumRecord& a, const SpectrumRecord& b)
    {
      return sameValue(a.target_mz, b.target_mz) && sameValue(a.lower_offset, b.lower_offset) &&
             sameValue(a.upper_offset, b.upper_offset) && sameValue(a.faims_cv, b.faims_cv) &&
             sameValue(a.ion_mobility_lower, b.ion_mobility_lower) && sameValue(a.ion_mobility_upper, b.ion_mobility_upper);
    };
    std::map<long long, std::vector<Size>> windows_by_target;
    std::vector<std::vector<const SpectrumRecord*>> window_records;
    std::vector<const SpectrumRecord*> without_window;
    for (const SpectrumRecord* r : ms2)
    {
      if (!r->has_isolation_window)
      {
        without_window.push_back(r);
        continue;
      }
      const long long bin = std::llround(r->target_mz * 1000.0);
      Size index = window_records.size();
      for (long long b = bin - 1; b <= bin + 1 && index == window_records.size(); ++b)
      {
        const auto candidates = windows_by_target.find(b);
        if (candidates == windows_by_target.end()) continue;
        for (Size w : candidates->second)
        {
          if (sameWindow(*window_records[w].front(), *r))
          {
            index = w;
            break;
          }
        }
      }
      if (index == window_records.size())
      {
        window_records.emplace_back();
        windows_by_target[bin].push_back(index);
      }
      window_records[index].push_back(r);
    }

    // the statistics of a group of MS2 spectra (sorted by retention time)
    auto groupMetrics = [](const std::vector<const SpectrumRecord*>& records, const GroupStatistics& stats)
    {
      WindowMetrics group;
      group.mass_resolving_power = stats.mass_resolving_power;
      group.spectrum_count = records.size();
      group.rt_min = stats.rt_min;
      group.rt_max = stats.rt_max;
      group.cycle_time_median = stats.cycle_time_median;
      group.tic_quantile_rt = stats.tic_quantile_rt;
      group.total_tic = stats.total_tic;
      group.peak_count = stats.peak_count;
      return group;
    };

    std::vector<double> pooled_rt_differences, window_cycle_times, window_spectra, window_lower, window_upper,
      window_widths, window_half_tic_rts, window_total_tics, window_peak_count_medians;
    for (const auto& records : window_records)
    {
      const SpectrumRecord& first = *records.front();
      const GroupStatistics stats(records);

      WindowMetrics window = groupMetrics(records, stats);
      window.has_isolation_window = true;
      window.target_mz = first.target_mz;
      window.lower_mz = first.target_mz - first.lower_offset;
      window.upper_mz = first.target_mz + first.upper_offset;
      window.width_mz = window.upper_mz - window.lower_mz;
      window.faims_cv = first.faims_cv;
      window.ion_mobility_lower = first.ion_mobility_lower;
      window.ion_mobility_upper = first.ion_mobility_upper;
      run.windows.push_back(window);

      if (window.spectrum_count == 1) ++run.windows_measured_once;
      pooled_rt_differences.insert(pooled_rt_differences.end(), stats.rt_differences.begin(), stats.rt_differences.end());
      window_cycle_times.push_back(window.cycle_time_median);
      window_spectra.push_back(static_cast<double>(window.spectrum_count));
      window_lower.push_back(window.lower_mz);
      window_upper.push_back(window.upper_mz);
      window_widths.push_back(window.width_mz);
      window_half_tic_rts.push_back(window.tic_quantile_rt[1]);
      window_total_tics.push_back(window.total_tic);
      window_peak_count_medians.push_back(window.peak_count[2]);
    }

    run.window_count = run.windows.size();
    if (!without_window.empty())
    {
      run.without_isolation_window = groupMetrics(without_window, GroupStatistics(without_window));
    }
    std::tie(run.window_spectra_min, run.window_spectra_max) = finiteRange(window_spectra);
    run.window_mz_min = finiteRange(window_lower).first;
    run.window_mz_max = finiteRange(window_upper).second;
    std::tie(run.window_width_min, run.window_width_max) = finiteRange(window_widths);
    std::tie(run.window_half_tic_rt_min, run.window_half_tic_rt_max) = finiteRange(window_half_tic_rts);
    std::tie(run.window_total_tic_min, run.window_total_tic_max) = finiteRange(window_total_tics);
    std::tie(run.window_peak_count_median_min, run.window_peak_count_median_max) = finiteRange(window_peak_count_medians);
    run.window_cycle_time_median = medianOrNaN(pooled_rt_differences);

    double sum = 0.0;
    Size n = 0;
    for (double t : window_cycle_times)
    {
      if (!std::isfinite(t)) continue;
      sum += t;
      ++n;
    }
    if (n > 0) run.window_cycle_time_mean = sum / static_cast<double>(n);

    return run;
  }

  void DIAQCMetrics::writeRunTable(const std::vector<RunMetrics>& runs, std::ostream& os)
  {
    writeRow(os, {"SourceFile", "Instrument", "SerialNumber", "StartTimeStamp", "RTDuration",
                  "mzMLMS1Count", "mzMLMSnCount", "MS1Resolution",
                  "MS1TIC25ileRT", "MS1TIC50ileRT", "MS1TIC75ileRT", "MS1TotalTIC", "MS1CycleTime",
                  "MS1PkCountMin", "MS1PkCount25ile", "MS1PkCount50ile", "MS1PkCount75ile", "MS1PkCountMax",
                  "IsolationWindowCount", "CyclesMin", "CyclesMax",
                  "MZRangeMin", "MZRangeMax", "IsolationWindowWidthMin", "IsolationWindowWidthMax",
                  "AverageMedianCycleTime", "TICMedianRTMin", "TICMedianRTMax", "TotalTICMin", "TotalTICMax",
                  "PkCountMedianMin", "PkCountMedianMax",
                  "MS2TIC25ileRT", "MS2TIC50ileRT", "MS2TIC75ileRT", "MS2TotalTIC",
                  "MS2PkCountMin", "MS2PkCount25ile", "MS2PkCount50ile", "MS2PkCount75ile", "MS2PkCountMax",
                  "MedianWindowCycleTime"});
    for (const RunMetrics& run : runs)
    {
      writeRow(os, {textCell(run.source_file), textCell(run.instrument), textCell(run.serial_number), textCell(run.start_time_stamp),
                    numberToString(toMinutes(run.rt_max)),
                    countToString(run.ms1_count), countToString(run.msn_count), numberToString(run.ms1_mass_resolving_power),
                    numberToString(toMinutes(run.ms1_tic_quantile_rt[0])), numberToString(toMinutes(run.ms1_tic_quantile_rt[1])),
                    numberToString(toMinutes(run.ms1_tic_quantile_rt[2])), numberToString(run.ms1_total_tic),
                    numberToString(run.ms1_cycle_time_median),
                    countToString(run.ms1_peak_count[0]), countToString(run.ms1_peak_count[1]), countToString(run.ms1_peak_count[2]),
                    countToString(run.ms1_peak_count[3]), countToString(run.ms1_peak_count[4]),
                    countToString(run.window_count), countToString(run.window_spectra_min), countToString(run.window_spectra_max),
                    numberToString(run.window_mz_min), numberToString(run.window_mz_max),
                    numberToString(run.window_width_min), numberToString(run.window_width_max),
                    numberToString(run.window_cycle_time_mean),
                    numberToString(toMinutes(run.window_half_tic_rt_min)), numberToString(toMinutes(run.window_half_tic_rt_max)),
                    numberToString(run.window_total_tic_min), numberToString(run.window_total_tic_max),
                    countToString(run.window_peak_count_median_min), countToString(run.window_peak_count_median_max),
                    numberToString(toMinutes(run.ms2_tic_quantile_rt[0])), numberToString(toMinutes(run.ms2_tic_quantile_rt[1])),
                    numberToString(toMinutes(run.ms2_tic_quantile_rt[2])), numberToString(run.ms2_total_tic),
                    countToString(run.ms2_peak_count[0]), countToString(run.ms2_peak_count[1]), countToString(run.ms2_peak_count[2]),
                    countToString(run.ms2_peak_count[3]), countToString(run.ms2_peak_count[4]),
                    numberToString(run.window_cycle_time_median)});
    }
  }

  void DIAQCMetrics::writeWindowTable(const std::vector<RunMetrics>& runs, std::ostream& os)
  {
    writeRow(os, {"SourceFile", "LoMZ", "HiMZ", "WidthMZ", "IonMobility", "MassResolvingPower",
                  "MSMSCount", "RTMin", "RTMax", "CycleTimeMedian", "TIC25ileRT", "TIC50ileRT", "TIC75ileRT",
                  "TotalTIC", "PkCountMin", "PkCount25ile", "PkCount50ile", "PkCount75ile", "PkCountMax",
                  "IonMobilityLow", "IonMobilityHigh", "TargetMZ"});
    for (const RunMetrics& run : runs)
    {
      std::vector<const WindowMetrics*> rows;
      for (const WindowMetrics& w : run.windows) rows.push_back(&w);
      if (run.without_isolation_window.spectrum_count > 0) rows.push_back(&run.without_isolation_window); // m/z values NA
      for (const WindowMetrics* row : rows)
      {
        const WindowMetrics& w = *row;
        writeRow(os, {textCell(run.source_file),
                      numberToString(w.lower_mz), numberToString(w.upper_mz), numberToString(w.width_mz),
                      numberToString(w.faims_cv), numberToString(w.mass_resolving_power),
                      countToString(w.spectrum_count), numberToString(toMinutes(w.rt_min)), numberToString(toMinutes(w.rt_max)),
                      numberToString(w.cycle_time_median),
                      numberToString(toMinutes(w.tic_quantile_rt[0])), numberToString(toMinutes(w.tic_quantile_rt[1])),
                      numberToString(toMinutes(w.tic_quantile_rt[2])), numberToString(w.total_tic),
                      countToString(w.peak_count[0]), countToString(w.peak_count[1]), countToString(w.peak_count[2]),
                      countToString(w.peak_count[3]), countToString(w.peak_count[4]),
                      numberToString(w.ion_mobility_lower), numberToString(w.ion_mobility_upper),
                      numberToString(w.target_mz)});
      }
    }
  }

  void DIAQCMetrics::writeMzQC(const std::vector<RunMetrics>& runs, std::ostream& os, const std::string& software_version, const std::string& creation_date)
  {
    using json = nlohmann::ordered_json;
    const ControlledVocabulary& cv = ControlledVocabulary::getPSIMSCV();
    std::set<std::string> reported_missing;

    // a CV term as {accession, name[, description]}; empty if the term is not in the vocabulary
    auto term = [&cv, &reported_missing](const std::string& accession, bool with_description = false) -> json
    {
      if (!cv.exists(accession))
      {
        if (reported_missing.insert(accession).second)
        {
          OPENMS_LOG_WARN << "Warning: CV term '" << accession << "' is not in the installed PSI-MS vocabulary. "
                          << "Values that need it are not written to the mzQC file." << std::endl;
        }
        return json();
      }
      const ControlledVocabulary::CVTerm& cv_term = cv.getTerm(accession);
      json result{{"accession", accession}, {"name", cv_term.name}};
      if (with_description && !cv_term.description.empty()) result["description"] = cv_term.description;
      return result;
    };

    auto allFinite = [](const std::vector<double>& values)
    {
      return std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); });
    };

    json run_qualities = json::array();
    for (const RunMetrics& run : runs)
    {
      json metrics = json::array();
      auto add = [&](const std::string& accession, const json& value, const std::string& unit)
      {
        json metric = term(accession, true);
        if (metric.is_null()) return;
        metric["value"] = value;
        if (!unit.empty())
        {
          json unit_term = term(unit);
          if (!unit_term.is_null()) metric["unit"] = unit_term;
        }
        metrics.push_back(metric);
      };
      // numbers that may be undefined: the metric is left out rather than written as null
      auto addNumbers = [&](const std::string& accession, const std::vector<double>& values, const std::string& unit, bool integer)
      {
        if (values.empty() || !allFinite(values)) return;
        json value = json::array();
        for (double v : values)
        {
          if (integer) value.push_back(std::llround(v));
          else value.push_back(v);
        }
        add(accession, values.size() == 1 ? value[0] : value, unit);
      };
      const std::string second = "UO:0000010", minute = "UO:0000031", count = "UO:0000189", mz = "MS:1000040", intensity = "MS:1000043";
      const auto& q1 = run.ms1_peak_count;
      const auto& q2 = run.ms2_peak_count;

      addNumbers("MS:4000067", {run.rt_max - run.rt_min}, second, false); // MS run duration
      addNumbers("MS:4000070", {run.rt_min, run.rt_max}, second, false); // retention time acquisition range
      add("MS:4000059", run.ms1_count, count); // number of MS1 spectra
      add("MS:4000060", run.ms2_count, count); // number of MS2 spectra
      addNumbers("MS:4000061", {q1[1], q1[2], q1[3]}, count, true); // MS1 density quantiles
      addNumbers("MS:4000062", {q2[1], q2[2], q2[3]}, count, true); // MS2 density quantiles
      addNumbers("MS:4000190", {toMinutes(run.ms1_tic_quantile_rt[0]), toMinutes(run.ms1_tic_quantile_rt[1]), toMinutes(run.ms1_tic_quantile_rt[2])}, minute, false); // MS1 TIC quantile RT
      addNumbers("MS:4000191", {toMinutes(run.ms2_tic_quantile_rt[0]), toMinutes(run.ms2_tic_quantile_rt[1]), toMinutes(run.ms2_tic_quantile_rt[2])}, minute, false); // MS2 TIC quantile RT
      addNumbers("MS:4000192", {run.ms1_cycle_time_median}, second, false); // MS1 median cycle time
      // m/z acquisition range: the range of the precursor m/z values of MSn spectra (not the m/z range that the
      // isolation windows cover, which is in the run table)
      addNumbers("MS:4000069", {run.precursor_mz_min, run.precursor_mz_max}, mz, false);
      if (run.window_count > 0) // only DIA isolation windows; MS2 spectra without an isolation window are not one
      {
        addNumbers("MS:4000193", {run.window_cycle_time_median}, second, false); // DIA isolation window median cycle time
        add("MS:4000194", run.window_count, count); // DIA isolation window count
        addNumbers("MS:4000195", {run.window_width_min, run.window_width_max}, mz, false); // DIA isolation window m/z widths
        // times a window is measured; note that PSI-MS 4.2.2 names this term like MS:4000194 ("DIA isolation window count")
        addNumbers("MS:4000196", {run.window_spectra_min, run.window_spectra_max}, count, true);
        addNumbers("MS:4000197", {toMinutes(run.window_half_tic_rt_min), toMinutes(run.window_half_tic_rt_max)}, minute, false); // DIA isolation window half TIC RT
        addNumbers("MS:4000198", {run.window_total_tic_min, run.window_total_tic_max}, intensity, false); // DIA isolation window TIC
        addNumbers("MS:4000199", {run.window_peak_count_median_min, run.window_peak_count_median_max}, count, true); // DIA isolation window peak count
      }

      json input_file;
      input_file["location"] = fileURI(run.input_path);
      input_file["name"] = File::basename(run.input_path);
      const auto [format_accession, format_name] = fileFormatTerm(run.input_type);
      input_file["fileFormat"] = json{{"accession", format_accession}, {"name", format_name}};
      json properties = json::array();
      auto addProperty = [&](const std::string& accession, const std::string& value)
      {
        if (value.empty()) return;
        json property = term(accession);
        if (property.is_null()) return;
        property["value"] = value;
        properties.push_back(property);
      };
      addProperty("MS:1000569", run.file_sha1); // SHA-1
      addProperty("MS:1000031", run.instrument); // instrument model
      addProperty("MS:1000529", run.serial_number); // instrument serial number
      if (!properties.empty()) input_file["fileProperties"] = properties;

      json software = term("MS:1000752", true); // TOPP software
      if (software.is_null()) software = json{{"accession", "MS:1000752"}, {"name", "TOPP software"}};
      software["value"] = "DIAuditor";
      software["version"] = software_version;
      software["uri"] = "https://www.openms.de";

      json run_quality;
      run_quality["metadata"]["label"] = run.source_file;
      run_quality["metadata"]["inputFiles"] = json::array({input_file});
      run_quality["metadata"]["analysisSoftware"] = json::array({software});
      run_quality["qualityMetrics"] = metrics;
      run_qualities.push_back(run_quality);
    }

    // PSI-MS carries the unit (UO) terms it uses; declaring UO as well would make them ambiguous for validators.
    // The URI names the release of the installed vocabulary (the purl always points to the latest one). The terms are
    // looked up in getPSIMSCV(), which the mzML reader has loaded already; its version() is that of the last OBO file
    // it loaded (not psi-ms.obo), so the version is read from the header of psi-ms.obo.
    const std::string ms_version = oboDataVersion(File::find("/CV/psi-ms.obo"));
    json vocabulary_ms{{"name", "Proteomics Standards Initiative Mass Spectrometry Ontology"}};
    if (ms_version.empty())
    {
      vocabulary_ms["uri"] = "http://purl.obolibrary.org/obo/ms/psi-ms.obo";
    }
    else
    {
      vocabulary_ms["uri"] = "https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v" + ms_version + "/psi-ms.obo";
      vocabulary_ms["version"] = ms_version;
    }

    json out;
    out["mzQC"]["version"] = "1.0.0";
    out["mzQC"]["creationDate"] = creation_date;
    out["mzQC"]["description"] = "Data-independent acquisition quality metrics (DIAuditor)";
    out["mzQC"]["runQualities"] = run_qualities;
    out["mzQC"]["controlledVocabularies"] = json::array({vocabulary_ms});
    // replace invalid UTF-8 (e.g. from a file name in another encoding) instead of failing
    os << out.dump(2, ' ', false, json::error_handler_t::replace) << '\n';
  }
} // namespace OpenMS
