// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <limits>
#include <set>
#include <string>

#include <openms_thermo_bridge/thermo_bridge.hpp>

namespace OpenMS
{
  namespace
  {
    /// Case-insensitive "haystack contains needle".
    bool containsCI(const std::string& haystack, const std::string& needle)
    {
      if (needle.empty()) { return true; }
      auto it = std::search(haystack.begin(), haystack.end(),
                            needle.begin(), needle.end(),
                            [](char a, char b)
                            { return std::tolower(static_cast<unsigned char>(a)) ==
                                     std::tolower(static_cast<unsigned char>(b)); });
      return it != haystack.end();
    }

    /// Parse a double, returning false if the string is empty/unparseable.
    bool parseDouble(const std::string& s, double& out)
    {
      if (s.empty()) { return false; }
      try
      {
        size_t pos = 0;
        out = std::stod(s, &pos);
        // Reject trailing garbage (e.g. "12.3abc") but tolerate surrounding whitespace.
        while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) { ++pos; }
        return pos == s.size();
      }
      catch (...)
      {
        return false;
      }
    }

    /// Find the value of the first trailer-extra label that contains @p needle (case-insensitive).
    /// The bridge matches keys exactly, so we pass the exact label string back to it.
    std::string trailerValueContaining(const openms::thermo_bridge::RawFile& raw, int scan,
                                       const std::vector<std::string>& labels, const std::string& needle)
    {
      for (const std::string& label : labels)
      {
        if (containsCI(label, needle))
        {
          return raw.trailer_extra_value(scan, label);
        }
      }
      return "";
    }

    /// Map a Thermo MassAnalyzerType enum string (e.g. "MassAnalyzerFTMS") to an OpenMS analyzer type.
    /// FTMS is disambiguated into Orbitrap vs. FT-ICR via the instrument model, mirroring
    /// the logic in openms-thermo-bridge cv_mapping::cv_mass_analyzer().
    MassAnalyzer::AnalyzerType mapAnalyzer(const std::string& t, const std::string& model)
    {
      if (t == "MassAnalyzerFTMS" || t == "FTMS")
      {
        if (containsCI(model, "Orbitrap") || containsCI(model, "Exactive") ||
            containsCI(model, "Exploris") || containsCI(model, "Astral"))
        {
          return MassAnalyzer::AnalyzerType::ORBITRAP;
        }
        return MassAnalyzer::AnalyzerType::FOURIERTRANSFORM; // FT-ICR
      }
      if (t == "MassAnalyzerITMS") { return MassAnalyzer::AnalyzerType::IT; }            // ion trap (MS:1000264)
      if (t == "MassAnalyzerTQMS" || t == "MassAnalyzerSQMS") { return MassAnalyzer::AnalyzerType::QUADRUPOLE; }
      if (t == "MassAnalyzerTOFMS" || t == "MassAnalyzerASTMS") { return MassAnalyzer::AnalyzerType::TOF; }
      if (t == "MassAnalyzerSector") { return MassAnalyzer::AnalyzerType::SECTOR; }
      return MassAnalyzer::AnalyzerType::ANALYZERNULL;
    }

    /// Map a Thermo IonizationModeType enum string (e.g. "ElectroSpray") to an OpenMS ionization method.
    IonSource::IonizationMethod mapIonization(const std::string& t)
    {
      if (t == "ElectroSpray") { return IonSource::IonizationMethod::ESI; }
      if (t == "NanoSpray" || t == "CardNanoSprayIonization") { return IonSource::IonizationMethod::NESI; }
      if (t == "AtmosphericPressureChemicalIonization") { return IonSource::IonizationMethod::APCI; }
      if (t == "ChemicalIonization") { return IonSource::IonizationMethod::CI; }
      if (t == "MatrixAssistedLaserDesorptionIonization") { return IonSource::IonizationMethod::MALDI; }
      if (t == "ElectronImpact") { return IonSource::IonizationMethod::EI; }
      if (t == "FastAtomBombardment") { return IonSource::IonizationMethod::FAB; }
      if (t == "ThermoSpray") { return IonSource::IonizationMethod::TSP; }
      if (t == "FieldDesorption") { return IonSource::IonizationMethod::FD; }
      if (t == "GlowDischarge") { return IonSource::IonizationMethod::GD_MS; }
      return IonSource::IonizationMethod::IONMETHODNULL;
    }

    /// Infer the ion detector type from the mass analyzer: Orbitrap/FT-ICR use image-current
    /// (inductive) detection, while ion traps and quadrupoles use an electron multiplier.
    /// (The bridge's cv_detector_types() is not exported, so we derive this locally.)
    IonDetector::Type detectorForAnalyzer(MassAnalyzer::AnalyzerType analyzer)
    {
      switch (analyzer)
      {
        case MassAnalyzer::AnalyzerType::ORBITRAP:
        case MassAnalyzer::AnalyzerType::FOURIERTRANSFORM:
          return IonDetector::Type::INDUCTIVEDETECTOR;
        case MassAnalyzer::AnalyzerType::IT:
        case MassAnalyzer::AnalyzerType::QUADRUPOLE:
          return IonDetector::Type::ELECTRONMULTIPLIER;
        default:
          return IonDetector::Type::TYPENULL;
      }
    }

    /// Extract the scan window (m/z range) from a Thermo scan filter string, e.g.
    /// "FTMS + p NSI Full ms [375.0000-1500.0000]" -> (375.0, 1500.0). Returns false if no
    /// single bracketed range is found (e.g. multi-range SIM filters with a comma are skipped).
    bool parseScanWindowFromFilter(const std::string& filter, double& lower, double& upper)
    {
      const std::string::size_type close = filter.rfind(']');
      if (close == std::string::npos) { return false; }
      const std::string::size_type open = filter.rfind('[', close);
      if (open == std::string::npos) { return false; }

      const std::string inside = filter.substr(open + 1, close - open - 1);
      if (inside.find(',') != std::string::npos) { return false; } // multiple ranges (e.g. SIM)
      const std::string::size_type dash = inside.find('-');
      if (dash == std::string::npos) { return false; }

      if (!parseDouble(inside.substr(0, dash), lower)) { return false; }
      if (!parseDouble(inside.substr(dash + 1), upper)) { return false; }
      return upper > lower && lower > 0.0;
    }
  } // anonymous namespace

  void ThermoRawFile::load(const std::string& path, MSExperiment& exp)
  {
    // Loaders are expected to populate the output experiment from scratch.
    exp = MSExperiment();

    if (!File::exists(path))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);
    }

    std::filesystem::path raw_path(static_cast<std::string>(path));

    try
    {
      openms::thermo_bridge::RawFile raw(raw_path);

      // --- Source file metadata ---
      SourceFile src;
      src.setNameOfFile(raw_path.filename().string());
      src.setPathToFile(raw_path.parent_path().string());
      src.setFileType("Thermo RAW format");
      src.setNativeIDType("Thermo nativeID format");
      src.setNativeIDTypeAccession("MS:1000768"); // Thermo nativeID format
      exp.getSourceFiles().push_back(src);

      // --- Instrument / experiment-level metadata ---
      const std::string model = raw.instrument_model();
      if (!model.empty())
      {
        exp.getInstrument().setName(raw.instrument_name());
        exp.getInstrument().setModel(model);
        exp.getInstrument().setCustomizations(raw.instrument_serial_number());
      }
      const std::string sw_version = raw.instrument_software_version();
      if (!sw_version.empty())
      {
        exp.getInstrument().getSoftware().setVersion(sw_version);
      }
      const std::string sample_name = raw.sample_name();
      if (!sample_name.empty())
      {
        exp.getSample().setName(sample_name);
      }
      const std::string creation_date = raw.creation_date(); // ISO-8601 (e.g. "2018-05-31T12:34:56.789...")
      if (creation_date.size() >= 19)
      {
        // OpenMS DateTime only parses up to seconds; drop fractional seconds / timezone offset.
        try
        {
          exp.setDateTime(DateTime::fromString(creation_date.substr(0, 19), "yyyy-MM-ddThh:mm:ss"));
        }
        catch (const Exception::BaseException&)
        {
          // leave date-time unset if the format is unexpected
        }
      }
      const std::string file_description = raw.file_description();
      if (!file_description.empty())
      {
        exp.getInstrument().setMetaValue("file description", file_description);
      }

      // --- Read spectra ---
      const int scan_count = raw.scan_count();
      const int first_scan = raw.first_scan_number();
      const int last_scan = raw.last_scan_number();

      OPENMS_LOG_INFO << "ThermoRawFile: reading " << scan_count << " scans from "
                      << raw_path.filename().string() << "\n";

      setProgress(0);
      startProgress(first_scan, last_scan, "Loading Thermo RAW file");

      exp.reserve(scan_count);

      // Distinct mass-analyzer strings and the first ionization mode seen, used to build the
      // single OpenMS Instrument component lists after the scan loop.
      std::set<std::string> analyzer_types;
      std::string ionization_mode;

      for (int scan = first_scan; scan <= last_scan; ++scan)
      {
        setProgress(scan);

        MSSpectrum spectrum;

        // RT in seconds (bridge returns minutes)
        const double rt_minutes = raw.retention_time(scan);
        spectrum.setRT(rt_minutes * 60.0);

        // MS level
        const int ms_level = raw.ms_level(scan);
        spectrum.setMSLevel(ms_level);

        // Native ID
        spectrum.setNativeID("scan=" + std::to_string(scan));

        // Centroid / profile
        const bool centroid = raw.is_centroid_scan(scan);
        spectrum.setType(centroid ? SpectrumSettings::SpectrumType::CENTROID : SpectrumSettings::SpectrumType::PROFILE);

        // Polarity — bridge returns the raw Thermo PolarityType enum value. Verified empirically
        // against real RAW data: Negative = 0, Positive = 1, Any = 2 (the bridge's own header comment
        // and cv_polarity() helper claim the opposite and are wrong).
        const int polarity = raw.polarity(scan);
        if (polarity == 1)
        {
          spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE);
        }
        else if (polarity == 0)
        {
          spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE);
        }

        // Scan filter string (PSI-MS MS:1000512). Also collect the analyzer type and the
        // configured scan window (m/z range) carried inside the filter.
        const std::string filter = raw.scan_filter(scan);
        if (!filter.empty())
        {
          spectrum.setMetaValue("filter string", filter);

          double win_lower = 0.0;
          double win_upper = 0.0;
          if (parseScanWindowFromFilter(filter, win_lower, win_upper))
          {
            ScanWindow window;
            window.begin = win_lower;
            window.end = win_upper;
            spectrum.getInstrumentSettings().getScanWindows().push_back(window);
          }
        }

        const std::string analyzer = raw.mass_analyzer_type(scan);
        if (!analyzer.empty()) { analyzer_types.insert(analyzer); }
        if (ionization_mode.empty()) { ionization_mode = raw.ionization_mode(scan); }

        // Scan statistics (PSI-MS): total ion current, base peak m/z and intensity.
        spectrum.setMetaValue("total ion current", raw.tic(scan));
        spectrum.setMetaValue("base peak m/z", raw.base_peak_mass(scan));
        spectrum.setMetaValue("base peak intensity", raw.base_peak_intensity(scan));

        // Trailer-extra metadata: ion injection time and FAIMS compensation voltage. Neither has a
        // dedicated bridge accessor; both are discovered by matching the (raw) trailer labels.
        const std::vector<std::string> trailer_labels = raw.trailer_extra_labels(scan);
        if (!trailer_labels.empty())
        {
          double injection_time = 0.0;
          if (parseDouble(trailerValueContaining(raw, scan, trailer_labels, "Ion Injection Time"), injection_time))
          {
            // Ion injection time (MS:1000927) is a <scan>-level term; storing it on an Acquisition
            // makes it serialize as a scan cvParam (matching ProteoWizard/msconvert output).
            Acquisition acq;
            acq.setMetaValue("ion injection time", injection_time);
            spectrum.getAcquisitionInfo().push_back(acq);
          }

          double faims_cv = 0.0;
          if (parseDouble(trailerValueContaining(raw, scan, trailer_labels, "FAIMS CV"), faims_cv))
          {
            spectrum.setDriftTime(faims_cv);
            spectrum.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
          }
        }

        // Precursor info for MSn
        if (ms_level > 1)
        {
          Precursor prec;
          const double prec_mz = raw.precursor_mass(scan);
          if (prec_mz > 0.0)
          {
            prec.setMZ(prec_mz);
          }
          const int charge = raw.precursor_charge(scan);
          if (charge != 0)
          {
            prec.setCharge(charge);
          }
          // Map activation type string to OpenMS activation method (independent of collision energy,
          // so electron-based methods that report no eV still get an activation method).
          const std::string act_type = raw.activation_type(scan);
          if (act_type == "CID" || act_type == "CollisionInducedDissociation")
          {
            prec.getActivationMethods().insert(Precursor::ActivationMethod::CID);
          }
          else if (act_type == "HCD" || act_type == "HigherEnergyCollisionalDissociation")
          {
            prec.getActivationMethods().insert(Precursor::ActivationMethod::HCD);
          }
          else if (act_type == "ETD" || act_type == "ElectronTransferDissociation")
          {
            prec.getActivationMethods().insert(Precursor::ActivationMethod::ETD);
          }
          else if (act_type == "ECD" || act_type == "ElectronCaptureDissociation")
          {
            prec.getActivationMethods().insert(Precursor::ActivationMethod::ECD);
          }
          else if (act_type == "PQD")
          {
            prec.getActivationMethods().insert(Precursor::ActivationMethod::PQD);
          }
          else if (!act_type.empty())
          {
            OPENMS_LOG_WARN << "ThermoRawFile: unknown activation type '" << act_type
                            << "' for scan " << scan << ", defaulting to CID\n";
            prec.getActivationMethods().insert(Precursor::ActivationMethod::CID);
          }
          const double ce = raw.collision_energy(scan);
          if (ce > 0.0)
          {
            prec.setActivationEnergy(ce);
          }
          const double iso_width = raw.isolation_width(scan);
          if (iso_width > 0.0)
          {
            prec.setIsolationWindowLowerOffset(iso_width / 2.0);
            prec.setIsolationWindowUpperOffset(iso_width / 2.0);
          }
          spectrum.getPrecursors().push_back(prec);
        }

        // Read spectrum data (prefer centroided data for centroid scans)
        auto data = raw.spectrum_data(scan, centroid);

        // Populate peaks, tracking the lowest/highest observed m/z (PSI-MS MS:1000528/MS:1000527).
        spectrum.resize(data.mz.size());
        double lowest_mz = std::numeric_limits<double>::max();
        double highest_mz = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < data.mz.size(); ++i)
        {
          const double mz = data.mz[i];
          spectrum[i].setMZ(mz);
          spectrum[i].setIntensity(static_cast<Peak1D::IntensityType>(data.intensities[i]));
          lowest_mz = std::min(lowest_mz, mz);
          highest_mz = std::max(highest_mz, mz);
        }
        if (!data.mz.empty())
        {
          spectrum.setMetaValue("lowest observed m/z", lowest_mz);
          spectrum.setMetaValue("highest observed m/z", highest_mz);
        }

        exp.addSpectrum(std::move(spectrum));
      }

      endProgress();

      // --- Instrument component lists (single OpenMS Instrument) ---
      // OpenMS represents one instrument configuration per experiment, whereas Thermo files may use
      // several (e.g. FTMS + ion trap). We capture the union of analyzers seen plus the ion source and
      // detector(s) implied by the instrument model.
      Instrument& instrument = exp.getInstrument();
      Int component_order = 1;
      if (!ionization_mode.empty())
      {
        const IonSource::IonizationMethod method = mapIonization(ionization_mode);
        if (method != IonSource::IonizationMethod::IONMETHODNULL)
        {
          IonSource source;
          source.setIonizationMethod(method);
          source.setOrder(component_order++);
          instrument.getIonSources().push_back(source);
        }
      }
      std::set<IonDetector::Type> detector_types;
      for (const std::string& analyzer : analyzer_types)
      {
        const MassAnalyzer::AnalyzerType type = mapAnalyzer(analyzer, model);
        if (type == MassAnalyzer::AnalyzerType::ANALYZERNULL) { continue; }
        MassAnalyzer mass_analyzer;
        mass_analyzer.setType(type);
        mass_analyzer.setOrder(component_order++);
        instrument.getMassAnalyzers().push_back(mass_analyzer);

        const IonDetector::Type detector = detectorForAnalyzer(type);
        if (detector != IonDetector::Type::TYPENULL) { detector_types.insert(detector); }
      }
      for (const IonDetector::Type type : detector_types)
      {
        IonDetector ion_detector;
        ion_detector.setType(type);
        ion_detector.setOrder(component_order++);
        instrument.getIonDetectors().push_back(ion_detector);
      }

      // --- TIC chromatogram ---
      const openms::thermo_bridge::ChromatogramData tic_data = raw.chromatogram_data();
      if (!tic_data.times.empty())
      {
        MSChromatogram tic;
        tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
        tic.setNativeID("TIC");
        for (size_t i = 0; i < tic_data.times.size(); ++i)
        {
          tic.push_back(ChromatogramPeak(tic_data.times[i] * 60.0, tic_data.intensities[i]));
        }
        exp.addChromatogram(std::move(tic));
      }

      // Sort by RT
      exp.sortSpectra(true);
      exp.updateRanges();

      OPENMS_LOG_INFO << "ThermoRawFile: loaded " << exp.size() << " spectra\n";
    }
    catch (const openms::thermo_bridge::bridge_error& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        path, std::string("Thermo bridge error: ") + e.what());
    }
  }

} // namespace OpenMS

#endif // WITH_THERMO_RAW
