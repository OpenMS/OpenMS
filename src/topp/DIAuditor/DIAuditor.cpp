// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, David L. Tabb $
// --------------------------------------------------------------------------

#include "DIAQCMetrics.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <map>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DIAuditor DIAuditor

@brief Quality metrics for data-independent acquisition (DIA) runs, per run and per isolation window.

This tool re-implements <a href="https://github.com/dtabb73/DIAuditor">DIAuditor</a> by David L. Tabb. It reads
one or more mzML files and reports how each run was acquired: the MS1 sampling, the DIA isolation windows (count, m/z
range, widths, ion mobility, how often each window was measured and how fast) and how ion current and peak counts are
distributed over the run and over the windows. This helps to develop DIA methods and to recover the method of DIA
data acquired elsewhere.

<B>Inputs</B>: mzML files are read spectrum by spectrum, so their size does not matter. Depending on the build, Thermo
.raw files (the OpenMS Thermo reader, which needs a .NET 8 runtime) and Bruker timsTOF .d directories or .d.zip
archives (diaPASEF frames are split into their isolation windows) can be given directly. These readers load a whole run
into memory, so large runs need correspondingly much memory; converting them to mzML first avoids that. Spectra are
used as they are stored, profile or centroided, which matters for the peak counts. Values a reader does not provide are
NA: the mass resolving power (MS1Resolution, MassResolvingPower) for .raw and .d input, and the instrument, its serial
number and the start time for .d input.

MS2 spectra are grouped into isolation windows by the isolation window of their precursor (target m/z, lower and upper
offset) and by their ion mobility settings (see @p ion_mobility): the FAIMS compensation voltage, and the ion mobility
range of the window as it is written for diaPASEF frames ('ion mobility lower/upper limit', e.g. by
<tt>msconvert --combineIonMobilitySpectra</tt> or the OpenMS timsTOF reader). These values are matched with an absolute
tolerance of 1e-6 (as in OpenSWATH). Windows are listed in the order in which they are first acquired. MS2 spectra
without an isolation window (no precursor, or no isolation window offsets) are not a DIA isolation window: they are
listed together in the last row of a run in @p out_windows, without m/z values, and left out of the window count and
of all other window metrics. Spectra of MS level 3 and higher, and spectra without retention time, are counted but not
assigned to windows.

Peak count quartiles are taken as in DIAuditor: the values at positions n/4, n/2 and n/4 + n/2 (integer division) of
the n sorted counts.

The tool warns about data it cannot describe well: MS2 spectra that are single ion mobility scans (e.g. diaPASEF
converted without combining the scans of a frame; every scan then counts as a measurement of its window), MS2 spectra
with several precursors (multiplexed DIA, e.g. MSX; only the first isolation window is used), MS2 spectra without an
isolation window, and runs without isolation windows or whose windows are mostly measured once (e.g. DDA).

<B>Outputs</B> (at least one is required)
- @p out: one row per run, with the columns of DIAuditor's <tt>DIAuditor-byRun.tsv</tt>, followed by run-level
  MS2 statistics.
- @p out_windows: one row per isolation window of each run, with the columns of DIAuditor's
  <tt>DIAuditor-byIsolationWindow.tsv</tt> (IonMobility is the FAIMS compensation voltage), followed by the ion
  mobility range of the window and the isolation target m/z. MS2 spectra without an isolation window follow in one
  row with NA m/z values.
- @p out_mzqc: an mzQC file with one runQuality per run. It contains only metrics defined in the PSI-MS vocabulary
  (e.g. MS:4000190 to MS:4000199 for DIA), with the units the vocabulary defines. Per-window values are not
  part of it, because the vocabulary has no term for them. The DIA window metrics are written only for runs with
  isolation windows. MS:4000069 'm/z acquisition range' is the range of the precursor m/z values of the MSn spectra
  (as the vocabulary defines it); the run table's MZRangeMin/MZRangeMax is the m/z range the isolation windows cover.

In the tables, retention times are given in minutes, cycle times in seconds; values that are undefined (e.g. the
cycle time of a window measured only once) are written as NA.

<B>Differences to the original DIAuditor</B>
- Input files are given explicitly (@p in) instead of all mzML files in the current directory, and may also be Thermo
  .raw files or Bruker .d directories.
- Scan start times are converted from any unit the mzML reader supports; DIAuditor treats times written in seconds as
  minutes.
- A FAIMS compensation voltage is attributed to its own spectrum; DIAuditor can attribute the value of an MS1 spectrum to
  the MS2 spectrum before it.
- Windows are also separated by their isolation offsets and by the ion mobility range of diaPASEF frames (unless
  @p ion_mobility is set to 'faims' or 'none'); DIAuditor uses the isolation target and the FAIMS voltage only.
- Only MS2 spectra form windows; DIAuditor also puts MS3 and higher spectra into them.
- MS2 spectra without an isolation window are not a window; DIAuditor puts them into a window at m/z 0, which also
  extends its m/z range down to 0.
- Spectra without a 'total ion current' value get the sum of their intensities (see @p tic), and the instrument
  model is recognised from the whole PSI-MS vocabulary.
- Values that define a window are matched with a tolerance of 1e-6; DIAuditor requires exact equality.
- Medians of times (cycle times) are the usual median, i.e. the mean of the two middle values for an even count;
  DIAuditor takes the value at position n/2 of the n - 1 sorted time differences. MS1Resolution and a window's
  MassResolvingPower are the median over the respective spectra; DIAuditor takes the last MS1 spectrum and the first
  spectrum of the window.
- The sum of intensities (see @p tic) is accumulated in double precision.
- RTDuration is the retention time of the last spectrum of any MS level; DIAuditor uses the last MSn spectrum.
- The column header 'IsolationWidowWidthMax' of DIAuditor is spelled 'IsolationWindowWidthMax'.
- Runs and windows with very few spectra do not stop the tool; undefined values are reported as NA.
- mzQC values follow the definitions and units of the vocabulary: e.g. the run duration is reported as
  MS:4000067 'MS run duration' (last minus first scan, in seconds), and each metric appears once per run.
- Runs are labelled by their file name, which therefore has to be unique among the inputs.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_DIAuditor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_DIAuditor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDIAuditor :
  public TOPPBase
{
public:
  TOPPDIAuditor() :
    TOPPBase("DIAuditor", "Computes quality metrics of data-independent acquisition (DIA) runs, per run and per isolation window.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", ListUtils::create<std::string>(""),
                           "Input files, one per run: mzML, Thermo .raw and Bruker timsTOF .d (directory or .d.zip), "
                           "as supported by this build");
    StringList in_formats = {"mzML"};
#ifdef WITH_THERMO_RAW
    in_formats.push_back("raw");
#endif
#ifdef WITH_OPENTIMS
    in_formats.push_back("d");
#endif
    setValidFormats_("in", in_formats);
    registerOutputFile_("out", "<file>", "", "Table with one row per run", false);
    setValidFormats_("out", {"tsv"});
    registerOutputFile_("out_windows", "<file>", "", "Table with one row per isolation window of each run", false);
    setValidFormats_("out_windows", {"tsv"});
    registerOutputFile_("out_mzqc", "<file>", "", "mzQC file with the run metrics", false);
    setValidFormats_("out_mzqc", {"mzQC"});

    registerStringOption_("tic", "<choice>", "auto",
      "Total ion current of a spectrum. 'auto': the value in the file if present, otherwise the sum of intensities; "
      "'file': only the value in the file, spectra without it count with 0 (like the original DIAuditor); "
      "'computed': always the sum of intensities.", false);
    setValidStrings_("tic", {"auto", "file", "computed"});
    registerStringOption_("peak_count", "<choice>", "all",
      "Peaks counted per spectrum. 'all': all data points (like the original DIAuditor); 'nonzero': data points with "
      "an intensity above zero (differs for profile data).", false);
    setValidStrings_("peak_count", {"all", "nonzero"});
    registerStringOption_("ion_mobility", "<choice>", "auto",
      "Ion mobility settings that separate isolation windows with the same m/z range. 'auto': the FAIMS compensation "
      "voltage and the ion mobility range of the window (diaPASEF frames); 'faims': only the FAIMS compensation voltage "
      "(like the original DIAuditor); 'none': windows are defined by m/z alone.", false);
    setValidStrings_("ion_mobility", {"auto", "faims", "none"});
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in = getStringList_("in");
    const std::string out = getStringOption_("out");
    const std::string out_windows = getStringOption_("out_windows");
    const std::string out_mzqc = getStringOption_("out_mzqc");
    if (out.empty() && out_windows.empty() && out_mzqc.empty())
    {
      OPENMS_LOG_ERROR << "Error: no output given. Set at least one of 'out', 'out_windows' and 'out_mzqc'." << endl;
      return ILLEGAL_PARAMETERS;
    }

    DIAQCMetrics::Options options;
    const std::string tic = getStringOption_("tic");
    options.tic_source = (tic == "file") ? DIAQCMetrics::TICSource::FILE :
                         (tic == "computed") ? DIAQCMetrics::TICSource::COMPUTED : DIAQCMetrics::TICSource::AUTO;
    options.peak_count = (getStringOption_("peak_count") == "nonzero") ? DIAQCMetrics::PeakCountMode::NONZERO :
                                                                       DIAQCMetrics::PeakCountMode::ALL;
    const std::string ion_mobility = getStringOption_("ion_mobility");
    options.ion_mobility = (ion_mobility == "faims") ? DIAQCMetrics::IonMobilityKey::FAIMS :
                           (ion_mobility == "none") ? DIAQCMetrics::IonMobilityKey::NONE : DIAQCMetrics::IonMobilityKey::AUTO;

    // paths without trailing directory separators (e.g. 'run.d/' from shell completion), as FileHandler::getType() does
    StringList files;
    for (std::string file : in)
    {
      while (file.size() > 1 && (file.ends_with('/') || file.ends_with('\\'))) file.pop_back();
      files.push_back(file);
    }

    // runs are labelled by file name, which has to be unique (also required by mzQC)
    std::map<std::string, std::string> labels;
    for (const std::string& file : files)
    {
      const auto [it, inserted] = labels.emplace(File::stemName(file), file);
      if (!inserted)
      {
        OPENMS_LOG_ERROR << "Error: '" << it->second << "' and '" << file << "' have the same file name. "
                         << "Runs are labelled by file name, so input file names must be unique." << endl;
        return ILLEGAL_PARAMETERS;
      }
    }

    std::vector<DIAQCMetrics::RunMetrics> runs;
    for (const std::string& file : files)
    {
      OPENMS_LOG_INFO << "Reading " << file << " ..." << endl;
      const FileTypes::Type type = FileHandler::getType(file);
      DIAQCMetrics metrics(options);
      if (type == FileTypes::MZML)
      {
        // stream the file: only a small record per spectrum is kept
        MSDataTransformingConsumer consumer;
        consumer.setExperimentalSettingsFunc([&metrics](const ExperimentalSettings& settings) { metrics.setExperimentalSettings(settings); });
        consumer.setSpectraProcessingFunc([&metrics](MSSpectrum& spectrum) { metrics.addSpectrum(spectrum); });
        MzMLFile mzml;
        mzml.setLogType(log_type_);
        mzml.getOptions().setSortSpectraByMZ(false); // peaks are only counted and summed
        mzml.transform(file, &consumer, true);
      }
      else
      {
        // the vendor readers (Thermo .raw, Bruker .d) load the whole run
        PeakMap exp;
        FileHandler().loadExperiment(file, exp, {type}, log_type_);
        metrics.setExperimentalSettings(exp);
        for (const MSSpectrum& spectrum : exp) metrics.addSpectrum(spectrum);
      }

      DIAQCMetrics::RunMetrics run = metrics.compute();
      run.source_file = File::stemName(file);
      run.input_path = file;
      run.input_type = type;
      // a .d directory has no single file to checksum
      if (!out_mzqc.empty() && !File::isDirectory(file)) run.file_sha1 = FileHandler::computeFileHash(file);

      OPENMS_LOG_INFO << File::basename(file) << ": " << run.ms1_count << " MS1 and " << run.msn_count
                      << " MSn spectra, " << run.window_count << " isolation windows." << endl;
      if (run.ms2_count == 0)
      {
        OPENMS_LOG_WARN << "Warning: no MS2 spectra in '" << file << "', so there are no isolation windows. Is this a DIA run?" << endl;
      }
      else if (run.window_count == 0 && run.without_isolation_window.spectrum_count > 0) // (MS2 spectra without RT: warned below)
      {
        OPENMS_LOG_WARN << "Warning: no MS2 spectrum of '" << file << "' has an isolation window, so there are no DIA "
                        << "isolation windows. Is this a DIA run?" << endl;
      }
      else if (run.windows_measured_once > run.window_count / 2)
      {
        OPENMS_LOG_WARN << "Warning: most isolation windows of '" << file << "' are measured only once. Is this a DIA run?" << endl;
      }
      if (run.window_count > 0 && run.without_isolation_window.spectrum_count > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << run.without_isolation_window.spectrum_count << " MS2 spectra of '" << file
                        << "' have no isolation window. They are not counted as a DIA isolation window." << endl;
      }
      if (run.spectra_without_rt > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << run.spectra_without_rt << " spectra of '" << file << "' have no retention time. "
                        << "They are counted, but not used for any other metric." << endl;
      }
      if (run.ms2_multiple_precursors > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << run.ms2_multiple_precursors << " MS2 spectra of '" << file << "' have more than "
                        << "one precursor (multiplexed DIA?). Only the first isolation window of each is used." << endl;
      }
      if (run.ms2_scan_ion_mobility > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << run.ms2_scan_ion_mobility << " MS2 spectra of '" << file << "' have an ion mobility "
                        << "of their own, i.e. they are single ion mobility scans. Each of them counts as a measurement of its "
                        << "window, so counts and cycle times describe scans, not frames. Combine the scans of a frame "
                        << "(e.g. msconvert --combineIonMobilitySpectra) for per-frame metrics." << endl;
      }
      runs.push_back(std::move(run));
    }

    auto open = [](const std::string& filename, std::ofstream& os)
    {
      os.open(filename);
      if (!os) throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    };
    if (!out.empty())
    {
      std::ofstream os;
      open(out, os);
      DIAQCMetrics::writeRunTable(runs, os);
    }
    if (!out_windows.empty())
    {
      std::ofstream os;
      open(out_windows, os);
      DIAQCMetrics::writeWindowTable(runs, os);
    }
    if (!out_mzqc.empty())
    {
      // mzQC requires an RFC 3339 date-time, i.e. with time zone
      std::ostringstream mzqc;
      DIAQCMetrics::writeMzQC(runs, mzqc, VersionInfo::getVersion(), DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ"));
      std::ofstream os;
      open(out_mzqc, os);
      os << mzqc.str();
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDIAuditor tool;
  return tool.main(argc, argv);
}

/// @endcond
