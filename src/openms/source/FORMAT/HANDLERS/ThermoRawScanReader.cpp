// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ThermoRawScanReader.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/FileNameUtils.h>
#include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <exception>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <nlohmann/json.hpp>
#include <openms_thermo_bridge/cv_mapping.hpp>
#include <openms_thermo_bridge/thermo_bridge.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS::Internal
{
namespace
{
  using Json = nlohmann::json;

  /// Case-insensitive "haystack contains needle".
  bool containsCI(const std::string& haystack, const std::string& needle)
  {
    if (needle.empty()) { return true; }
    auto it = std::search(haystack.begin(), haystack.end(), needle.begin(), needle.end(),
                          [](char a, char b) { return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b)); });
    return it != haystack.end();
  }

  /// True if @p dir holds the two files the bridge needs to start the .NET runtime.
  bool isManagedDirectory(const std::filesystem::path& dir)
  {
    std::error_code ec;
    return std::filesystem::is_regular_file(dir / "ThermoWrapperManaged.dll", ec)
        && std::filesystem::is_regular_file(dir / "ThermoWrapperManaged.runtimeconfig.json", ec);
  }

  /**
    @brief Locate the managed half of the thermo bridge (ThermoWrapperManaged.dll + runtimeconfig).

    The bridge's own default looks next to its shared library (<lib dir>/managed or
    <lib dir>/openms_thermo_bridge/managed). That works for a plain OpenMS build or install
    tree, but not for relocated layouts such as Python wheels, where wheel repair tools
    rename and move the bridge library. OpenMS therefore checks, in order:

      1. the OPENMS_THERMO_MANAGED_DIR environment variable, if it holds the managed files
         (a stale or mistyped override is reported and skipped rather than breaking a
         correctly bundled installation),
      2. <OpenMS share dir>/openms_thermo_bridge/managed, i.e. the copy that OpenMS installs
         into its shared-data directory (pyOpenMS ships share/OpenMS inside the wheel),
      3. nothing, so the bridge falls back to its own lookup.
  */
  std::optional<std::filesystem::path> resolveManagedDirectory()
  {
    if (const char* env = std::getenv("OPENMS_THERMO_MANAGED_DIR"); env != nullptr && env[0] != '\0')
    {
      std::filesystem::path dir(env);
      if (isManagedDirectory(dir))
      {
        return dir;
      }
      OPENMS_LOG_WARN << "ThermoRawFile: OPENMS_THERMO_MANAGED_DIR='" << env
                      << "' does not contain ThermoWrapperManaged.dll and its runtimeconfig.json; "
                      << "ignoring it and looking in the OpenMS share directory instead\n";
    }

    std::filesystem::path share_dir(File::getOpenMSDataPath());
    std::filesystem::path candidate = share_dir / "openms_thermo_bridge" / "managed";
    if (isManagedDirectory(candidate))
    {
      return candidate;
    }
    return std::nullopt;
  }

  /// Map a Thermo MassAnalyzerType enum string (e.g. "MassAnalyzerFTMS") to an
  /// OpenMS analyzer type. FTMS is disambiguated into Orbitrap vs. FT-ICR via the
  /// instrument model, mirroring the logic in openms-thermo-bridge
  /// cv_mapping::cv_mass_analyzer().
  MassAnalyzer::AnalyzerType mapAnalyzer(const std::string& t, const std::string& model)
  {
    if (t == "MassAnalyzerFTMS" || t == "FTMS")
    {
      if (containsCI(model, "Orbitrap") || containsCI(model, "Exactive") || containsCI(model, "Exploris") || containsCI(model, "Astral"))
      {
        return MassAnalyzer::AnalyzerType::ORBITRAP;
      }
      return MassAnalyzer::AnalyzerType::FOURIERTRANSFORM; // FT-ICR
    }
    if (t == "MassAnalyzerITMS") { return MassAnalyzer::AnalyzerType::IT; } // ion trap (MS:1000264)
    if (t == "MassAnalyzerTQMS" || t == "MassAnalyzerSQMS") { return MassAnalyzer::AnalyzerType::QUADRUPOLE; }
    if (t == "MassAnalyzerTOFMS" || t == "MassAnalyzerASTMS") { return MassAnalyzer::AnalyzerType::TOF; }
    if (t == "MassAnalyzerSector") { return MassAnalyzer::AnalyzerType::SECTOR; }
    return MassAnalyzer::AnalyzerType::ANALYZERNULL;
  }

  /// Map a Thermo IonizationModeType enum string (e.g. "ElectroSpray") to an
  /// OpenMS ionization method.
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

  /// Infer the ion detector type from the mass analyzer: Orbitrap/FT-ICR use
  /// image-current (inductive) detection, while ion traps and quadrupoles use an
  /// electron multiplier. Choose the detector for this scan's analyzer, rather
  /// than every detector in the instrument.
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

  std::string text(const Json& object, const char* key)
  {
    const auto it = object.find(key);
    return it != object.end() && it->is_string() ? it->get<std::string>() : std::string();
  }

  double numeric(const Json& object, const char* key, double fallback = 0.0)
  {
    const auto it = object.find(key);
    return it != object.end() && it->is_number() ? it->get<double>() : fallback;
  }

  bool isNumber(const Json& object, const char* key)
  {
    const auto it = object.find(key);
    return it != object.end() && it->is_number();
  }

  DataValue unitValue(double value, int unit)
  {
    DataValue result(value);
    result.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
    result.setUnit(unit);
    return result;
  }

  DataValue msUnitValue(double value, int unit)
  {
    DataValue result(value);
    result.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
    result.setUnit(unit);
    return result;
  }

  void copyValues(const Json& object, MetaInfoInterface& target, const std::string& prefix)
  {
    for (const auto& [key, value] : object.items())
    {
      if (value.is_number_integer()) { target.setMetaValue(prefix + key, value.get<Int64>()); }
      else if (value.is_number()) { target.setMetaValue(prefix + key, value.get<double>()); }
      else if (value.is_string()) { target.setMetaValue(prefix + key, value.get<std::string>()); }
    }
  }

  namespace tb = openms::thermo_bridge;

  /// The vendor value, or unset if it is not finite (the bridge keeps NaN for missing values).
  std::optional<double> finite(double value)
  {
    return std::isfinite(value) ? std::optional<double>(value) : std::nullopt;
  }

  /**
    Typed view of one scan header for the precursor parser. Trailer values are the
    vendor-formatted strings (ReadOptions::trailer_text) paired with the labels of the
    file's trailer layout; reaction values that are not finite stay unset.
  */
  void toThermoScan(const tb::ScanHeader& header, const tb::TrailerLayout& layout, ThermoScan& scan)
  {
    scan.scan_number = header.scan_number;
    scan.ms_level = header.ms_level();
    scan.filter = header.filter;
    scan.native_id = header.native_id();
    const Size entries = std::min(layout.labels.size(), header.trailer_text.size());
    scan.trailer.resize(entries);
    for (Size i = 0; i < entries; ++i)
    {
      scan.trailer[i].label = layout.labels[i].label;
      scan.trailer[i].value = header.trailer_text[i];
    }
    scan.reactions.clear();
    for (const auto& reaction : header.reactions)
    {
      ThermoReaction converted;
      converted.precursor_mass = finite(reaction.precursor_mass);
      converted.isolation_width = finite(reaction.isolation_width);
      converted.isolation_offset = finite(reaction.isolation_offset);
      converted.activation = tb::activation_type_name(reaction.activation);
      converted.collision_energy = finite(reaction.collision_energy);
      converted.collision_energy_valid = reaction.collision_energy_valid;
      scan.reactions.push_back(std::move(converted));
    }
  }

  /// Append @p value as a JSON string literal, exactly as nlohmann::json::dump() writes it.
  void appendJsonString(std::string& out, const std::string& value)
  {
    const bool plain = std::all_of(value.begin(), value.end(), [](char c) {
      return c >= 0x20 && c < 0x7F && c != '"' && c != '\\';
    });
    if (plain)
    {
      out += '"';
      out += value;
      out += '"';
    }
    else
    {
      out += Json(value).dump(); // escapes and validates UTF-8
    }
  }

  /**
    The trailer as the JSON array of {label, value} objects that 'Thermo trailer extra' holds;
    same text as dumping the equivalent nlohmann::json array, without building it per scan.
  */
  std::string trailerJson(const ThermoScan& scan)
  {
    std::string out = "[";
    for (const auto& entry : scan.trailer)
    {
      if (out.size() > 1) { out += ','; }
      out += "{\"label\":";
      appendJsonString(out, entry.label);
      out += ",\"value\":";
      appendJsonString(out, entry.value);
      out += '}';
    }
    out += ']';
    return out;
  }

  /// Annotate @p precursor with the activation of @p reaction.
  void setActivation(Precursor& precursor, const ThermoReaction& reaction, bool supplemental)
  {
    const std::string& type = reaction.activation;
    if (reaction.collision_energy_valid && reaction.collision_energy)
    {
      precursor.setMetaValue(supplemental ? "supplemental collision energy" : "collision energy", unitValue(*reaction.collision_energy, 266));
    }
    if (supplemental)
    {
      // Keep the vendor terms and record the combined method (EThcD / ETciD) that downstream
      // consumers such as TheoreticalSpectrumGenerator use; MzMLHandler writes the same pair.
      if (type == "HigherEnergyCollisionalDissociation")
      {
        precursor.setMetaValue("supplemental beam-type collision-induced dissociation", "");
        precursor.getActivationMethods().insert(Precursor::ActivationMethod::EThcD);
      }
      else if (type == "CollisionInducedDissociation")
      {
        precursor.setMetaValue("supplemental collision-induced dissociation", "");
        precursor.getActivationMethods().insert(Precursor::ActivationMethod::ETciD);
      }
      else
      {
        precursor.setMetaValue("Thermo supplemental activation", type);
      }
      return;
    }
    static const std::map<std::string, Precursor::ActivationMethod> methods = {
      {"CollisionInducedDissociation", Precursor::ActivationMethod::CID}, {"HigherEnergyCollisionalDissociation", Precursor::ActivationMethod::HCD},
      {"ElectronTransferDissociation", Precursor::ActivationMethod::ETD}, {"ElectronCaptureDissociation", Precursor::ActivationMethod::ECD},
      {"MultiPhotonDissociation", Precursor::ActivationMethod::IMD},      {"PQD", Precursor::ActivationMethod::PQD}};
    auto it = methods.find(type);
    if (it != methods.end()) { precursor.getActivationMethods().insert(it->second); }
    else
    {
      const auto term = openms::thermo_bridge::cv_activation_type(type);
      if (term.accession != "MS:1000044") { precursor.setMetaValue(term.name, ""); }
      precursor.setMetaValue("Thermo activation type", type);
    }
  }

  /// Copy bridge peaks into @p spectrum (replacing its peaks) and sort them by m/z if needed.
  void fillPeaks(MSSpectrum& spectrum, const std::vector<double>& mz, const std::vector<double>& intensities)
  {
    spectrum.resize(mz.size());
    for (Size i = 0; i < mz.size(); ++i)
    {
      spectrum[i].setMZ(mz[i]);
      spectrum[i].setIntensity(intensities[i]);
    }
    if (! spectrum.isSorted()) { spectrum.sortByPosition(); }
  }

  /// Per-thread state of the scan conversion (the trailer index caches label positions).
  struct ScanConverter
  {
    ThermoRawScanReader::Options options;
    DataProcessingPtr processing;
    ThermoTrailerIndex trailer;
    ThermoTrailerIndex::Key injection_key = trailer.add("Ion Injection Time (ms):");
    ThermoTrailerIndex::Key mono_key = trailer.add("Monoisotopic M/Z:");
    ThermoTrailerIndex::Key orbitrap_resolution_key = trailer.add("Orbitrap Resolution:");
    ThermoTrailerIndex::Key ft_resolution_key = trailer.add("FT Resolution:");
    ThermoTrailerIndex::Key faims_on_key = trailer.add("FAIMS Voltage On:");
    ThermoTrailerIndex::Key faims_cv_key = trailer.add("FAIMS CV:");
  };

  /**
    @brief Convert one scan; shared by the sequential reader and the random access accessors.

    @param[in] c Conversion state
    @param[in] scan Scan number
    @param[in] h Scan header as read with trailer text
    @param[in] parsed Typed view of @p h
    @param[in] precursors Precursors resolved in acquisition order (empty for neutral loss/gain scans)
    @param[in] configuration_id Instrument configuration of the scan
    @param[in] first_seen First scan whose precursors were resolved (parents before it are unknown)
    @param[in] parent Peaks of a parent scan for the precursor intensity estimate
    @param[in] data Peaks and auxiliary arrays (nullptr: metadata only, no intensity estimate)
    @param[out] spectrum The converted spectrum
  */
  void convertScan(ScanConverter& c, int scan, const tb::ScanHeader& h, const ThermoScan& parsed, const std::vector<ThermoPrecursor>& precursors,
                   const std::string& configuration_id, int first_seen, const std::function<const MSSpectrum&(int)>& parent,
                   const tb::Spectrum* data, MSSpectrum& spectrum)
  {
    const int level = std::max(h.ms_level(), 0);
    const std::string order = tb::ms_order_name(h.ms_order);
    const bool neutral = h.ms_order == tb::MsOrder::Nl || h.ms_order == tb::MsOrder::Ng;
    spectrum = MSSpectrum();
    spectrum.setMSLevel(level);
    auto mode = level == 1 ? InstrumentSettings::ScanMode::MS1SPECTRUM : InstrumentSettings::ScanMode::MSNSPECTRUM;
    if (order == "Par") { mode = InstrumentSettings::ScanMode::PRECURSOR; }
    else if (order == "Nl") { mode = InstrumentSettings::ScanMode::CNL; }
    else if (order == "Ng") { mode = InstrumentSettings::ScanMode::CNG; }
    spectrum.getInstrumentSettings().setScanMode(mode);
    spectrum.setRT(h.retention_time * 60.0);
    spectrum.setNativeID(parsed.native_id);
    const bool centroid = c.options.centroid || h.centroid;
    spectrum.setType(centroid ? SpectrumSettings::SpectrumType::CENTROID : SpectrumSettings::SpectrumType::PROFILE);
    const int polarity = static_cast<int>(h.polarity);
    if (polarity == 1) { spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE); }
    else if (polarity == 0) { spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE); }
    Acquisition acquisition;
    acquisition.setMetaValue("filter string", parsed.filter);
    // Keep the historical spectrum-level accessor as well as the mzML scan-level term.
    spectrum.setMetaValue("filter string", parsed.filter);
    if (c.options.preserve_trailers) { acquisition.setMetaValue("Thermo trailer extra", trailerJson(parsed)); }
    const auto injection = ThermoRawFileMetadata::number(c.trailer.value(parsed, c.injection_key));
    if (injection) { acquisition.setMetaValue("ion injection time", unitValue(*injection, 28)); }
    const auto mono = ThermoRawFileMetadata::number(c.trailer.value(parsed, c.mono_key));
    if (mono && *mono > 0) { acquisition.setMetaValue("[Thermo Trailer Extra]Monoisotopic M/Z:", *mono); }
    // mass resolving power as msconvert takes it: 'Orbitrap Resolution:', else 'FT Resolution:' (older instruments).
    // Kept with the spectrum, as the mzML reader does; the mzML writer puts it into the scan (MS:1000800).
    auto resolution = ThermoRawFileMetadata::number(c.trailer.value(parsed, c.orbitrap_resolution_key));
    if (! resolution || *resolution <= 0) { resolution = ThermoRawFileMetadata::number(c.trailer.value(parsed, c.ft_resolution_key)); }
    if (resolution && *resolution > 0) { spectrum.setMetaValue("mass resolving power", std::llround(*resolution)); }
    const std::string& enabled_text = c.trailer.value(parsed, c.faims_on_key);
    const auto voltage_on = ThermoRawFileMetadata::number(enabled_text);
    const auto voltage = ThermoRawFileMetadata::number(c.trailer.value(parsed, c.faims_cv_key));
    std::string enabled = enabled_text;
    std::transform(enabled.begin(), enabled.end(), enabled.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (((voltage_on && *voltage_on != 0) || (enabled == "true" || enabled == "on" || enabled == "yes")) && voltage)
    {
      spectrum.setDriftTime(*voltage);
      spectrum.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
    }
    if (configuration_id != "ic_0") { acquisition.setMetaValue("instrument_configuration_ref", configuration_id); }
    spectrum.getAcquisitionInfo().push_back(acquisition);
    ScanWindow window;
    window.begin = h.low_mass;
    window.end = h.high_mass;
    if (std::isfinite(h.low_mass) && std::isfinite(h.high_mass) && window.end >= window.begin)
    {
      spectrum.getInstrumentSettings().getScanWindows().push_back(window);
    }
    if (std::isfinite(h.tic)) { spectrum.setMetaValue("total ion current", h.tic); }
    if (std::isfinite(h.base_peak_mass)) { spectrum.setMetaValue("base peak m/z", msUnitValue(h.base_peak_mass, 1000040)); }
    if (std::isfinite(h.base_peak_intensity))
    {
      spectrum.setMetaValue("base peak intensity", msUnitValue(h.base_peak_intensity, 1000131));
    }
    if (neutral && ! parsed.reactions.empty())
    {
      const ThermoReaction& reaction = parsed.reactions.front();
      spectrum.setMetaValue(order == "Nl" ? "neutral loss" : "neutral gain", reaction.precursor_mass.value_or(0.0));
      Product product;
      const double width = reaction.isolation_width.value_or(-1.0);
      const double offset = reaction.isolation_offset.value_or(0.0);
      if (width >= 0 && width / 2 >= std::abs(offset))
      {
        product.setIsolationWindowLowerOffset(width / 2 - offset);
        product.setIsolationWindowUpperOffset(width / 2 + offset);
      }
      spectrum.getProducts().push_back(product);
    }
    for (const auto& descriptor : precursors)
    {
      Precursor precursor;
      precursor.setMZ(descriptor.selected_mz); // OpenMS default: selected ion
      precursor.setMetaValue("isolation window target m/z", descriptor.target_mz);
      precursor.setMetaValue("selected ion m/z", descriptor.selected_mz);
      if (descriptor.charge) { precursor.setCharge(*descriptor.charge); }
      if (descriptor.width)
      {
        const double lower = descriptor.lower_offset, upper = descriptor.upper_offset;
        if (lower >= 0) { precursor.setIsolationWindowLowerOffset(lower); }
        if (upper >= 0) { precursor.setIsolationWindowUpperOffset(upper); }
        if (lower <= 0) { precursor.setMetaValue("isolation window lower offset", lower); }
        if (upper <= 0) { precursor.setMetaValue("isolation window upper offset", upper); }
      }
      if (! descriptor.spectrum_ref.empty()) { precursor.setMetaValue("spectrum_ref", descriptor.spectrum_ref); }
      setActivation(precursor, descriptor.activation, false);
      if (descriptor.supplemental) { setActivation(precursor, *descriptor.supplemental, true); }
      // Parents are scans whose header was read before this one.
      const bool known_parent = descriptor.parent_scan >= first_seen && descriptor.parent_scan < scan;
      if (data != nullptr && descriptor.estimate_intensity && known_parent && precursor.getMZ() > 0)
      {
        // TRFP's precursor-intensity estimate sums the target +/- 1.5 m/z.
        double intensity = 0;
        const MSSpectrum& parent_spectrum = parent(descriptor.parent_scan);
        const double target = descriptor.target_mz;
        const double half_width = descriptor.width.value_or(0.0) > 0 ? 1.5 : 0.0;
        for (auto peak = parent_spectrum.MZBegin(target - half_width); peak != parent_spectrum.end() && peak->getMZ() < target + half_width; ++peak)
        {
          intensity += peak->getIntensity();
        }
        precursor.setIntensity(intensity);
        if (intensity == 0) { precursor.setMetaValue("peak intensity", 0.0); }
        precursor.setMetaValue("peak intensity unit accession", "MS:1000131");
      }
      spectrum.getPrecursors().push_back(precursor);
    }
    spectrum.getDataProcessing().push_back(c.processing);
    if (data == nullptr) { return; }
    if (! data->mz.empty())
    {
      if (centroid)
      {
        const Size basepeak = std::distance(data->intensities.begin(), std::max_element(data->intensities.begin(), data->intensities.end()));
        spectrum.setMetaValue("base peak m/z", msUnitValue(data->mz[basepeak], 1000040));
        spectrum.setMetaValue("base peak intensity", msUnitValue(data->intensities[basepeak], 1000131));
      }
      // The bridge returns peaks sorted by m/z.
      spectrum.setMetaValue("lowest observed m/z", msUnitValue(data->mz.front(), 1000040));
      spectrum.setMetaValue("highest observed m/z", msUnitValue(data->mz.back(), 1000040));
    }
    spectrum.resize(data->mz.size());
    for (Size i = 0; i < data->mz.size(); ++i)
    {
      spectrum[i].setMZ(data->mz[i]);
      spectrum[i].setIntensity(data->intensities[i]);
    }
    if (! data->charges.empty() && data->charges.size() == spectrum.size())
    {
      MSSpectrum::IntegerDataArray array;
      array.setName("charge array");
      array.reserve(data->charges.size());
      for (double charge : data->charges)
      {
        array.push_back(static_cast<Int64>(charge));
      }
      spectrum.getIntegerDataArrays().push_back(std::move(array));
    }
    // The independently sampled noise grid is not a peak annotation; store double lists,
    // which MzMLHandler serializes as 64-bit binary arrays.
    if (c.options.noise_data)
    {
      if (! data->noise_mz.empty()) { spectrum.setMetaValue("sampled noise m/z array", data->noise_mz); }
      if (! data->noise_intensities.empty()) { spectrum.setMetaValue("sampled noise intensity array", data->noise_intensities); }
      if (! data->noise_baselines.empty()) { spectrum.setMetaValue("sampled noise baseline array", data->noise_baselines); }
    }
  }

  /// Shared state of the random access accessors: the handle that new accessors are created from.
  struct AccessorRoot
  {
    AccessorRoot(tb::RawFile h, const ScanConverter& c, int first, std::shared_ptr<const TempDir> u) :
      unpacked(std::move(u)), handle(std::move(h)), converter(c), first_seen(first)
    {
    }
    std::shared_ptr<const TempDir> unpacked; ///< RAW file extracted from an archive; declared first, so it outlives handle
    std::mutex mutex; ///< guards handle
    tb::RawFile handle;
    ScanConverter converter;
    int first_seen = 0;
  };
} // namespace

struct ThermoRawScanReader::Impl
{
  /// Owns the RAW file extracted from a ZIP archive (null otherwise). Shared with the accessors, which may outlive the
  /// reader; declared first, so it is removed after every handle on the file has been closed.
  std::shared_ptr<const TempDir> unpacked;
  Options options;
  std::filesystem::path raw_path; ///< as given by the caller (the archive for an extracted file)
  std::optional<openms::thermo_bridge::RawFile> raw;
  Json metadata;
  std::string sha1; ///< of the RAW file, only with Options::checksum
  Instrument instrument;
  std::string model;
  DataProcessingPtr processing;
  std::map<std::string, std::string> configurations; ///< analyzer:ionization -> instrument configuration id
  std::vector<std::pair<std::string, Instrument>> new_configurations; ///< seen by header(), not yet in the settings
  std::string configuration_id; ///< configuration of the scan whose header was read last

  ThermoRawFileMetadata precursor_parser;
  ScanConverter converter; ///< options and data processing are set by the constructor

  // State of the scan whose header was read last.
  ScanInfo info;
  tb::Spectrum header_read; ///< header (trailer text, peaks if requested) of the last header() call
  bool header_has_peaks = false; ///< header_read also holds the peaks
  tb::Spectrum peak_read;   ///< peaks read by convert() (buffers reused)

  /// What convert() reads for the peaks and auxiliary arrays.
  tb::ReadOptions peakReadOptions() const
  {
    tb::ReadOptions read;
    read.centroid = options.centroid;
    read.trailers = false;
    read.charges = options.charge_data && options.centroid;
    read.noise = options.noise_data;
    return read;
  }
  /// Instrument configuration id of a scan; a new analyzer / ionization combination is registered.
  const std::string& configurationId(const tb::ScanHeader& h)
  {
    const std::string analyzer = tb::mass_analyzer_name(h.analyzer), ionization = tb::ionization_mode_name(h.ionization);
    const std::string configuration_key = analyzer + ":" + ionization;
    auto configuration = configurations.find(configuration_key);
    if (configuration == configurations.end())
    {
      const std::string id = configurations.empty() ? "ic_0" : "thermo_ic_" + std::to_string(configurations.size());
      configuration = configurations.emplace(configuration_key, id).first;
      Instrument config = instrument;
      IonSource ion_source;
      ion_source.setIonizationMethod(mapIonization(ionization));
      ion_source.setOrder(1);
      ion_source.setMetaValue("ionization accession", openms::thermo_bridge::cv_ionization_mode(ionization).accession);
      config.getIonSources().push_back(ion_source);
      MassAnalyzer mass_analyzer;
      mass_analyzer.setType(mapAnalyzer(analyzer, model));
      mass_analyzer.setOrder(2);
      mass_analyzer.setMetaValue("mass analyzer accession", openms::thermo_bridge::cv_mass_analyzer(analyzer, model).accession);
      config.getMassAnalyzers().push_back(mass_analyzer);
      IonDetector detector;
      detector.setType(analyzer == "MassAnalyzerASTMS" ? IonDetector::Type::CONVERSIONDYNODEELECTRONMULTIPLIER
                                                       : detectorForAnalyzer(mass_analyzer.getType()));
      detector.setOrder(3);
      config.getIonDetectors().push_back(detector);
      new_configurations.emplace_back(id, std::move(config));
    }
    return configuration->second;
  }

  ThermoScan scan;
  std::string order;
  bool neutral = false;
  std::vector<ThermoPrecursor> precursors;
  int first_seen = 0; ///< first scan whose header was read (all later ones up to info.scan were read too)

  /**
    @brief Reads scans ahead on worker threads, so vendor decoding overlaps with the conversion.

    Workers claim chunks of consecutive scans and read each scan with its own bridge accessor;
    the consumer takes them in scan order. At most @p window chunks are read ahead of the one
    being consumed. Errors are rethrown when the consumer reaches the failed chunk.
  */
  struct ReadAhead
  {
    static constexpr int chunk = 128;
    int first = 0, last = 0, peaks_from = 0;
    tb::ReadOptions header_options, full_options;
    std::vector<std::thread> workers;
    std::mutex mutex;
    std::condition_variable produced, consumed;
    std::map<int, std::pair<std::vector<tb::Spectrum>, std::exception_ptr>> ready; ///< chunk -> scans
    int next_chunk = 0, current_chunk = -1, window = 0;
    bool stop = false;
    std::vector<tb::Spectrum> current;

    void run(tb::RawFile accessor)
    {
      for (;;)
      {
        int c;
        {
          std::unique_lock<std::mutex> lock(mutex);
          consumed.wait(lock, [&] { return stop || next_chunk < std::max(current_chunk, 0) + window; });
          if (stop) { return; }
          c = next_chunk++;
        }
        const int begin = first + c * chunk;
        if (begin > last) { return; }
        const int end = std::min(last, begin + chunk - 1);
        std::vector<tb::Spectrum> scans(end - begin + 1);
        std::exception_ptr error;
        try
        {
          for (int scan = begin; scan <= end; ++scan)
          {
            accessor.read_scan(scan, scans[scan - begin], scan >= peaks_from ? full_options : header_options);
          }
        }
        catch (...) { error = std::current_exception(); }
        {
          std::lock_guard<std::mutex> lock(mutex);
          ready.emplace(c, std::make_pair(std::move(scans), error));
        }
        produced.notify_all();
      }
    }

    /// Swap the read-ahead data of @p scan into @p out; false if @p scan is not covered.
    bool take(int scan, tb::Spectrum& out)
    {
      if (scan < first || scan > last) { return false; }
      const int c = (scan - first) / chunk;
      if (c != current_chunk)
      {
        std::exception_ptr error;
        {
          std::unique_lock<std::mutex> lock(mutex);
          produced.wait(lock, [&] { return ready.count(c) > 0; });
          auto it = ready.find(c);
          current = std::move(it->second.first);
          error = it->second.second;
          ready.erase(it);
          current_chunk = c;
        }
        consumed.notify_all();
        if (error) { std::rethrow_exception(error); }
      }
      std::swap(out, current[scan - first - c * chunk]);
      return true;
    }

    void shutdown()
    {
      {
        std::lock_guard<std::mutex> lock(mutex);
        stop = true;
      }
      consumed.notify_all();
      for (auto& worker : workers) { worker.join(); }
      workers.clear();
    }

    ~ReadAhead() { shutdown(); }
  };
  std::unique_ptr<ReadAhead> read_ahead;

  /// Shared by the random access accessors (created by accessor()).
  mutable std::shared_ptr<AccessorRoot> accessor_root;
  mutable std::mutex accessor_mutex;

  /// Parent spectra decoded on demand (when the caller does not have them), most recent last.
  std::deque<std::pair<int, MSSpectrum>> parent_cache;

  const MSSpectrum* parentSpectrum(int parent_scan, const ParentLookup& lookup)
  {
    if (lookup)
    {
      if (const MSSpectrum* found = lookup(parent_scan)) { return found; }
    }
    for (const auto& [number, spectrum] : parent_cache)
    {
      if (number == parent_scan) { return &spectrum; }
    }
    auto data = raw->spectrum_data(parent_scan, options.centroid);
    MSSpectrum spectrum;
    fillPeaks(spectrum, data.mz, data.intensities);
    if (parent_cache.size() >= 8) { parent_cache.pop_front(); }
    parent_cache.emplace_back(parent_scan, std::move(spectrum));
    return &parent_cache.back().second;
  }
};

ThermoRawScanReader::ThermoRawScanReader(const std::string& path, const Options& options) :
  impl_(std::make_unique<Impl>())
{
  Impl& d = *impl_;
  d.options = options;
  if (! File::exists(path)) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path); }
  d.raw_path = std::filesystem::path(path);
  std::filesystem::path open_path;
  // Thermo's RawFileReader opens plain RAW files only. A ZIP archive that holds one .raw file is extracted into a
  // temporary directory, which is removed with the reader and its accessors (as BrukerTimsFile does for a .d.zip).
  // Other compression (gzip, bzip2) is rejected with a clear message instead of the vendor library's failure to open
  // the file. The source file metadata keeps the path given by the user.
  const FileTypes::Type compression = FileNameUtils::compressionType(path);
  if (compression == FileTypes::ZIP)
  {
    const std::string entry = ZipArchiveFile::findSingleEntry(path, ".raw");
    std::unique_ptr<TempDir> dir;
    open_path = std::filesystem::path(ZipArchiveFile::extractEntryToTempFile(path, entry, dir));
    d.unpacked = std::move(dir);
  }
  else if (compression != FileTypes::UNKNOWN)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path,
                                std::string(compression == FileTypes::GZ ? "A gzip" : "A bzip2")
                                + "-compressed Thermo RAW file cannot be read: decompress it first. (A ZIP archive "
                                  "that holds one .raw file can be read as it is.)");
  }
  else
  {
    // Thermo's RawFileReader cannot open a file through a symbolic link (as staged by Nextflow/Galaxy),
    // so hand it the resolved target.
    std::error_code canonical_error;
    open_path = std::filesystem::canonical(d.raw_path, canonical_error);
    if (canonical_error) { open_path = d.raw_path; }
  }
  const std::optional<std::filesystem::path> managed_dir = resolveManagedDirectory();
  if (managed_dir) { d.raw.emplace(open_path, *managed_dir); }
  else { d.raw.emplace(open_path); }
  d.metadata = Json::parse(d.raw->file_metadata_json(options.instrument_methods));
  if (options.checksum) { d.sha1 = d.raw->file_sha1(); }
  if (d.metadata.at("schema_version") != 1) { throw std::runtime_error("Unsupported Thermo metadata schema"); }

  if (d.metadata.at("instrument").is_object())
  {
    const auto& info = d.metadata.at("instrument");
    d.model = text(info, "model");
    d.instrument.setName(openms::thermo_bridge::cv_instrument_model(d.model).name);
    d.instrument.setModel(d.model);
    d.instrument.setMetaValue("Thermo instrument name", text(info, "name"));
    d.instrument.setMetaValue("instrument serial number", text(info, "serial_number"));
    d.instrument.setMetaValue("Thermo hardware version", text(info, "hardware_version"));
    d.instrument.setMetaValue("Thermo instrument units", text(info, "units"));
    d.instrument.getSoftware().setName("Thermo acquisition software");
    d.instrument.getSoftware().setVersion(text(info, "software_version"));
  }
  d.processing = std::make_shared<DataProcessing>();
  d.processing->getSoftware().setName("OpenMS Thermo RAW reader");
  d.processing->getSoftware().setVersion(VersionInfo::getVersion());
  d.processing->setMetaValue("Thermo RawFileReader version", text(d.metadata, "reader_version"));
  d.processing->getProcessingActions().insert(DataProcessing::ProcessingAction::FORMAT_CONVERSION);
  if (options.centroid) { d.processing->getProcessingActions().insert(DataProcessing::ProcessingAction::PEAK_PICKING); }
  d.converter.options = options;
  d.converter.processing = d.processing;
}

ThermoRawScanReader::~ThermoRawScanReader() = default;

bool ThermoRawScanReader::hasMSData() const
{
  return impl_->raw->has_ms_data();
}

int ThermoRawScanReader::firstScan() const
{
  return impl_->raw->first_scan_number();
}

int ThermoRawScanReader::lastScan() const
{
  return impl_->raw->last_scan_number();
}

void ThermoRawScanReader::fillSettings(ExperimentalSettings& settings) const
{
  const Impl& d = *impl_;
  const auto& file = d.metadata.at("file");
  SourceFile source;
  source.setNameOfFile(d.raw_path.filename().string());
  source.setPathToFile(d.raw_path.parent_path().string());
  source.setFileType("Thermo RAW format");
  source.setNativeIDType("Thermo nativeID format");
  source.setNativeIDTypeAccession("MS:1000768");
  if (! d.sha1.empty()) { source.setChecksum(d.sha1, SourceFile::ChecksumType::SHA1); }
  source.setMetaValue("RAW file revision", file.at("revision").get<int>());
  source.setMetaValue("file description", text(file, "description"));
  settings.getSourceFiles().push_back(source);
  const std::string date = text(file, "creation_date"); // ISO-8601, e.g. "2018-05-31T12:34:56.789Z"
  if (date.size() >= 19)
  {
    try
    {
      // OpenMS DateTime only resolves seconds; the full timestamp is kept for the mzML startTimeStamp.
      settings.setDateTime(DateTime::fromString(date.substr(0, 19), "yyyy-MM-ddThh:mm:ss"));
      settings.setMetaValue("mzml_start_time_stamp", date);
    }
    catch (const Exception::BaseException&)
    {
      OPENMS_LOG_WARN << "ThermoRawFile: could not parse creation date '" << date << "'\n";
    }
  }
  const auto& sample = d.metadata.at("sample");
  copyValues(sample, settings.getSample(), "Thermo ");
  settings.getSample().setName(text(sample, "sample name"));
  settings.getSample().setNumber(text(sample, "sample number"));
  settings.getSample().setComment(text(sample, "sample comment"));
  settings.getSample().setMetaValue("Thermo user sample fields", d.metadata.at("user_sample_fields").dump());
  if (d.options.instrument_methods) { settings.setMetaValue("Thermo instrument methods", d.metadata.at("instrument_methods").dump()); }
  if (d.metadata.at("run").is_object()) { copyValues(d.metadata.at("run"), settings, "Thermo run "); }
  settings.setMetaValue("Thermo RawFileReader version", text(d.metadata, "reader_version"));
  settings.setInstrument(d.instrument);
}

void ThermoRawScanReader::readAhead(int last_scan, int peaks_from)
{
  Impl& d = *impl_;
  if (d.read_ahead) { d.read_ahead->shutdown(); d.read_ahead.reset(); }
  int threads = 1;
#ifdef _OPENMP
  if (! omp_in_parallel()) { threads = omp_get_max_threads(); }
#endif
  // The vendor reader scales up to about three concurrent accessors; one thread stays for the conversion.
  const int workers = std::min(3, threads - 1);
  const int first = d.first_seen == 0 ? firstScan() : d.info.scan + 1;
  last_scan = std::min(last_scan, lastScan());
  if (workers < 1 || first > last_scan) { return; }
  auto ahead = std::make_unique<Impl::ReadAhead>();
  ahead->first = first;
  ahead->last = last_scan;
  ahead->peaks_from = peaks_from;
  ahead->window = 4 * workers;
  ahead->full_options = d.peakReadOptions();
  ahead->full_options.trailer_text = true;
  ahead->header_options = ahead->full_options;
  ahead->header_options.peaks = false;
  ahead->header_options.noise = false;
  std::vector<tb::RawFile> accessors;
  for (int i = 0; i < workers; ++i) { accessors.push_back(d.raw->thread_accessor()); }
  for (auto& accessor : accessors)
  {
    ahead->workers.emplace_back([state = ahead.get(), accessor = std::move(accessor)]() mutable { state->run(std::move(accessor)); });
  }
  d.read_ahead = std::move(ahead);
}

const ThermoRawScanReader::ScanInfo& ThermoRawScanReader::header(int scan, bool with_peaks)
{
  Impl& d = *impl_;
  if (d.first_seen == 0) { d.first_seen = scan; }
  else if (scan != d.info.scan + 1) { throw std::logic_error("ThermoRawScanReader::header() must visit scans in order"); }
  if (d.read_ahead && d.read_ahead->take(scan, d.header_read))
  {
    d.header_has_peaks = scan >= d.read_ahead->peaks_from;
  }
  else
  {
    tb::ReadOptions read = d.peakReadOptions();
    read.peaks = with_peaks;
    read.noise = with_peaks && read.noise;
    read.trailer_text = true;
    d.raw->read_scan(scan, d.header_read, read);
    d.header_has_peaks = with_peaks;
  }
  const tb::ScanHeader& h = d.header_read.header;
  toThermoScan(h, d.raw->trailer_layout(), d.scan);
  d.order = tb::ms_order_name(h.ms_order);
  d.neutral = h.ms_order == tb::MsOrder::Nl || h.ms_order == tb::MsOrder::Ng;
  // The precursor parser has to see every scan in acquisition order, also those that are not converted.
  if (d.neutral) { d.precursors.clear(); }
  else { d.precursors = d.precursor_parser.precursors(d.scan); }
  d.configuration_id = d.configurationId(h);
  const int level = h.ms_level();
  d.info.scan = scan;
  d.info.ms_level = level > 0 ? level : 0;
  d.info.rt = h.retention_time * 60.0;
  d.info.precursor_mz = d.precursors.empty() ? std::numeric_limits<double>::quiet_NaN() : d.precursors.front().selected_mz;
  d.info.tic = std::isfinite(h.tic) ? h.tic : std::numeric_limits<double>::quiet_NaN();
  return d.info;
}

void ThermoRawScanReader::scanHeaders(int first, int last, bool precursors, const std::function<bool(const ScanInfo&)>& visit)
{
  Impl& d = *impl_;
  tb::ReadOptions read;
  read.peaks = false;
  read.trailers = false;
  read.trailer_text = precursors;
  ThermoRawFileMetadata parser;
  ThermoScan parsed;
  tb::ScanBatch batch;
  ScanInfo info;
  // Chunks keep the batch buffers small while amortizing the per-call overhead of the bridge.
  constexpr int chunk = 256;
  for (int from = first; from <= last; from += chunk)
  {
    d.raw->read_scans(from, std::min(last, from + chunk - 1), batch, read);
    for (const tb::ScanHeader& h : batch.headers)
    {
      d.configurationId(h);
      const int level = h.ms_level();
      info.scan = h.scan_number;
      info.ms_level = level > 0 ? level : 0;
      info.rt = h.retention_time * 60.0;
      info.tic = std::isfinite(h.tic) ? h.tic : std::numeric_limits<double>::quiet_NaN();
      info.precursor_mz = std::numeric_limits<double>::quiet_NaN();
      if (precursors && h.ms_order != tb::MsOrder::Nl && h.ms_order != tb::MsOrder::Ng)
      {
        toThermoScan(h, d.raw->trailer_layout(), parsed);
        const auto resolved = parser.precursors(parsed);
        if (! resolved.empty()) { info.precursor_mz = resolved.front().selected_mz; }
      }
      if (! visit(info)) { return; }
    }
  }
}

void ThermoRawScanReader::addInstrumentConfigurations(ExperimentalSettings& settings)
{
  for (auto& [id, config] : impl_->new_configurations)
  {
    if (id == "ic_0") { settings.setInstrument(config); }
    else { settings.getInstrumentConfigurations()[id] = std::move(config); }
  }
  impl_->new_configurations.clear();
}

void ThermoRawScanReader::convert(MSSpectrum& spectrum, ExperimentalSettings& settings, const ParentLookup& parent, bool peaks)
{
  Impl& d = *impl_;
  // Ids follow the order of all scans read, so they do not depend on the scan selection.
  addInstrumentConfigurations(settings);
  const tb::Spectrum* data = nullptr;
  if (peaks)
  {
    data = &d.header_read;
    if (! d.header_has_peaks)
    {
      d.raw->read_scan(d.info.scan, d.peak_read, d.peakReadOptions());
      data = &d.peak_read;
    }
  }
  convertScan(d.converter, d.info.scan, d.header_read.header, d.scan, d.precursors, d.configuration_id, d.first_seen,
              [&](int parent_scan) -> const MSSpectrum& { return *d.parentSpectrum(parent_scan, parent); }, data, spectrum);
}

ThermoRawScanReader::ScanRecord ThermoRawScanReader::record() const
{
  const Impl& d = *impl_;
  return ScanRecord{d.info, d.precursors, d.configuration_id, d.scan.native_id};
}

struct ThermoRawScanReader::Accessor::Impl
{
  std::shared_ptr<AccessorRoot> root;
  tb::RawFile raw;
  ScanConverter converter;
  int first_seen = 0;
  tb::Spectrum read;    ///< buffers reused across scans
  tb::Spectrum parent_read;
  ThermoScan parsed;
  /// Recently decoded parent spectra, most recently used last.
  std::deque<std::pair<int, MSSpectrum>> parents;
  static constexpr Size parent_capacity = 8;

  explicit Impl(std::shared_ptr<AccessorRoot> r) :
    root(std::move(r)),
    raw([this] {
      std::lock_guard<std::mutex> lock(root->mutex);
      return root->handle.thread_accessor();
    }()),
    converter(root->converter),
    first_seen(root->first_seen)
  {
  }

  tb::ReadOptions peakReadOptions() const
  {
    tb::ReadOptions options;
    options.centroid = converter.options.centroid;
    options.trailers = false;
    options.charges = converter.options.charge_data && converter.options.centroid;
    options.noise = converter.options.noise_data;
    return options;
  }

  const MSSpectrum& parent(int scan)
  {
    for (auto it = parents.begin(); it != parents.end(); ++it)
    {
      if (it->first == scan)
      {
        if (std::next(it) != parents.end())
        {
          auto entry = std::move(*it);
          parents.erase(it);
          parents.push_back(std::move(entry));
        }
        return parents.back().second;
      }
    }
    // Same decode as the sequential reader's fallback for parents the caller does not have.
    auto data = raw.spectrum_data(scan, converter.options.centroid);
    MSSpectrum spectrum;
    fillPeaks(spectrum, data.mz, data.intensities);
    if (parents.size() >= parent_capacity) { parents.pop_front(); }
    parents.emplace_back(scan, std::move(spectrum));
    return parents.back().second;
  }
};

ThermoRawScanReader::Accessor::Accessor(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}

ThermoRawScanReader::Accessor::~Accessor() = default;

void ThermoRawScanReader::Accessor::convert(const ScanRecord& record, MSSpectrum& spectrum, bool peaks)
{
  Impl& d = *impl_;
  tb::ReadOptions read = d.peakReadOptions();
  read.peaks = peaks;
  read.noise = peaks && read.noise;
  read.trailer_text = true;
  d.raw.read_scan(record.info.scan, d.read, read);
  toThermoScan(d.read.header, d.raw.trailer_layout(), d.parsed);
  convertScan(d.converter, record.info.scan, d.read.header, d.parsed, record.precursors, record.configuration_id, d.first_seen,
              [&d](int parent_scan) -> const MSSpectrum& { return d.parent(parent_scan); }, peaks ? &d.read : nullptr, spectrum);
}

void ThermoRawScanReader::Accessor::readPeaks(int scan, MSSpectrum& spectrum)
{
  Impl& d = *impl_;
  tb::ReadOptions read = d.peakReadOptions();
  read.noise = false;
  d.raw.read_scan(scan, d.read, read);
  spectrum = MSSpectrum();
  const tb::Spectrum& data = d.read;
  spectrum.resize(data.mz.size());
  for (Size i = 0; i < data.mz.size(); ++i)
  {
    spectrum[i].setMZ(data.mz[i]);
    spectrum[i].setIntensity(data.intensities[i]);
  }
  if (! data.charges.empty() && data.charges.size() == spectrum.size())
  {
    MSSpectrum::IntegerDataArray array;
    array.setName("charge array");
    array.reserve(data.charges.size());
    for (double charge : data.charges) { array.push_back(static_cast<Int64>(charge)); }
    spectrum.getIntegerDataArrays().push_back(std::move(array));
  }
}

std::unique_ptr<ThermoRawScanReader::Accessor> ThermoRawScanReader::Accessor::clone() const
{
  return std::unique_ptr<Accessor>(new Accessor(std::make_unique<Impl>(impl_->root)));
}

std::unique_ptr<ThermoRawScanReader::Accessor> ThermoRawScanReader::accessor() const
{
  const Impl& d = *impl_;
  std::shared_ptr<AccessorRoot> root;
  {
    std::lock_guard<std::mutex> lock(d.accessor_mutex);
    if (! d.accessor_root)
    {
      d.accessor_root = std::make_shared<AccessorRoot>(d.raw->thread_accessor(), d.converter, d.first_seen, d.unpacked);
    }
    root = d.accessor_root;
  }
  return std::unique_ptr<Accessor>(new Accessor(std::make_unique<Accessor::Impl>(root)));
}

std::shared_ptr<DataProcessing> ThermoRawScanReader::dataProcessing() const
{
  return impl_->processing;
}

MSChromatogram ThermoRawScanReader::ticChromatogram() const
{
  const auto data = impl_->raw->chromatogram_data();
  MSChromatogram tic;
  tic.setNativeID("TIC");
  tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
  tic.reserve(data.times.size());
  for (Size i = 0; i < data.times.size(); ++i)
  {
    tic.push_back(ChromatogramPeak(data.times[i] * 60.0, data.intensities[i]));
  }
  tic.getDataProcessing().push_back(impl_->processing);
  return tic;
}

void ThermoRawScanReader::readDetectors(std::vector<MSChromatogram>& chromatograms, std::vector<MSSpectrum>& spectra)
{
  Impl& d = *impl_;
  const Json detectors = Json::parse(d.raw->detector_chromatograms_json());
  for (const auto& trace : detectors.at("chromatograms"))
  {
    MSChromatogram chromatogram;
    chromatogram.setNativeID(text(trace, "native_id"));
    chromatogram.setName(text(trace, "label"));
    const std::string label = text(trace, "label"), device = text(trace, "device");
    const bool absorption = device == "UV" || device == "Pda";
    chromatogram.setChromatogramType(absorption ? ChromatogramSettings::ChromatogramType::ABSORPTION_CHROMATOGRAM
                                                : ChromatogramSettings::ChromatogramType::MASS_CHROMATOGRAM);
    chromatogram.setMetaValue("Thermo detector units", text(trace, "units"));
    if (absorption) { chromatogram.setMetaValue("mzml intensity array", "absorption"); }
    else if (containsCI(label, "pressure"))
    {
      chromatogram.setMetaValue("chromatogram type accession", "MS:1003019");
      chromatogram.setMetaValue("mzml intensity array", "pressure");
    }
    else if (containsCI(label, "flow"))
    {
      chromatogram.setMetaValue("chromatogram type accession", "MS:1003020");
      chromatogram.setMetaValue("mzml intensity array", "flow");
    }
    else if (! containsCI(label, "current") && ! containsCI(label, "fid"))
    {
      chromatogram.setMetaValue("chromatogram type accession", "MS:1000626");
      chromatogram.setMetaValue("mzml intensity array", "nonstandard");
    }
    const auto times = trace.at("times").get<std::vector<double>>();
    const auto values = trace.at("intensities").get<std::vector<double>>();
    if (times.size() != values.size()) { throw std::runtime_error("Inconsistent Thermo detector array lengths"); }
    chromatogram.reserve(times.size());
    for (Size i = 0; i < times.size(); ++i)
    {
      chromatogram.push_back(ChromatogramPeak(times[i] * 60.0, values[i]));
    }
    chromatogram.getDataProcessing().push_back(d.processing);
    chromatograms.push_back(std::move(chromatogram));
  }
  for (const auto& controller : d.metadata.at("controllers"))
  {
    if (text(controller, "name") != "Pda") { continue; }
    // A separate accessor keeps the MS controller selected for header() and convert().
    openms::thermo_bridge::RawFile pda = d.raw->thread_accessor();
    pda.select_instrument(controller.at("type").get<int>(), controller.at("number").get<int>());
    for (int scan = pda.first_scan_number(); scan <= pda.last_scan_number(); ++scan)
    {
      const Json meta = Json::parse(pda.scan_metadata_json(scan));
      MSSpectrum spectrum;
      spectrum.setNativeID(text(meta, "native_id"));
      spectrum.setRT(meta.at("retention_time").get<double>() * 60.0);
      spectrum.setMSLevel(0);
      spectrum.getInstrumentSettings().setScanMode(InstrumentSettings::ScanMode::ABSORPTION);
      spectrum.setType(SpectrumSettings::SpectrumType::PROFILE);
      spectrum.setMetaValue("mzml coordinate array", "wavelength");
      spectrum.setMetaValue("mzml intensity array", "absorption");
      Acquisition acquisition;
      spectrum.getAcquisitionInfo().push_back(acquisition);
      ScanWindow window;
      window.begin = numeric(meta, "low_wavelength");
      window.end = numeric(meta, "high_wavelength");
      window.setMetaValue("unit_accession", "UO:0000018");
      if (isNumber(meta, "low_wavelength") && isNumber(meta, "high_wavelength"))
      {
        spectrum.getInstrumentSettings().getScanWindows().push_back(window);
      }
      const auto data = pda.spectrum_data(scan, false);
      if (! data.mz.empty())
      {
        const auto [lo, hi] = std::minmax_element(data.mz.begin(), data.mz.end());
        spectrum.setMetaValue("lowest observed wavelength", unitValue(*lo, 18));
        spectrum.setMetaValue("highest observed wavelength", unitValue(*hi, 18));
      }
      fillPeaks(spectrum, data.mz, data.intensities);
      spectrum.getDataProcessing().push_back(d.processing);
      spectra.push_back(std::move(spectrum));
    }
  }
}
} // namespace OpenMS::Internal
