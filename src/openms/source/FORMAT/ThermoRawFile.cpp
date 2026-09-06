// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU
// Berlin SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

  #include <OpenMS/CONCEPT/Exception.h>
  #include <OpenMS/CONCEPT/LogStream.h>
  #include <OpenMS/DATASTRUCTURES/DateTime.h>
  #include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>
  #include <OpenMS/FORMAT/ThermoRawFile.h>
  #include <OpenMS/KERNEL/MSChromatogram.h>
  #include <OpenMS/KERNEL/MSSpectrum.h>
  #include <OpenMS/METADATA/Acquisition.h>
  #include <OpenMS/METADATA/AcquisitionInfo.h>
  #include <OpenMS/METADATA/DataProcessing.h>
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
  #include <map>
  #include <openms_thermo_bridge/cv_mapping.hpp>
  #include <openms_thermo_bridge/thermo_bridge.hpp>
  #include <set>
  #include <string>

namespace OpenMS
{
namespace
{
  /// Case-insensitive "haystack contains needle".
  bool containsCI(const std::string& haystack, const std::string& needle)
  {
    if (needle.empty()) { return true; }
    auto it = std::search(haystack.begin(), haystack.end(), needle.begin(), needle.end(),
                          [](char a, char b) { return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b)); });
    return it != haystack.end();
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

} // anonymous namespace

void ThermoRawFile::load(const std::string& path, MSExperiment& exp)
{
  exp = MSExperiment();
  if (! File::exists(path)) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path); }
  const std::filesystem::path raw_path(path);
  using Json = nlohmann::json;
  using Metadata = Internal::ThermoRawFileMetadata;
  auto text = [](const Json& object, const std::string& key) -> std::string {
    return object.contains(key) && object[key].is_string() ? object[key].get<std::string>() : "";
  };
  auto numeric = [](const Json& object, const std::string& key, double fallback = 0.0) -> double {
    return object.contains(key) && object[key].is_number() ? object[key].get<double>() : fallback;
  };
  auto unit_value = [](double value, int unit) -> DataValue {
    DataValue result(value);
    result.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
    result.setUnit(unit);
    return result;
  };
  auto ms_unit_value = [](double value, int unit) -> DataValue {
    DataValue result(value);
    result.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
    result.setUnit(unit);
    return result;
  };
  auto copy_values = [](const Json& object, MetaInfoInterface& target, const std::string& prefix = "") {
    for (const auto& [key, value] : object.items())
    {
      if (value.is_number_integer()) { target.setMetaValue(prefix + key, value.get<Int64>()); }
      else if (value.is_number()) { target.setMetaValue(prefix + key, value.get<double>()); }
      else if (value.is_string()) { target.setMetaValue(prefix + key, value.get<std::string>()); }
    }
  };
  try
  {
    openms::thermo_bridge::RawFile raw(raw_path);
    const Json metadata = Json::parse(raw.file_metadata_json(options_.instrument_methods, options_.checksum));
    if (metadata.at("schema_version") != 1) { throw std::runtime_error("Unsupported Thermo metadata schema"); }
    const auto& file = metadata.at("file");
    SourceFile source;
    source.setNameOfFile(raw_path.filename().string());
    source.setPathToFile(raw_path.parent_path().string());
    source.setFileType("Thermo RAW format");
    source.setNativeIDType("Thermo nativeID format");
    source.setNativeIDTypeAccession("MS:1000768");
    if (! text(file, "sha1").empty()) { source.setChecksum(text(file, "sha1"), SourceFile::ChecksumType::SHA1); }
    source.setMetaValue("RAW file revision", file.at("revision").get<int>());
    source.setMetaValue("file description", text(file, "description"));
    exp.getSourceFiles().push_back(source);
    const std::string date = text(file, "creation_date");
    if (date.size() >= 19)
    {
      exp.setDateTime(DateTime::fromString(date.substr(0, 19), "yyyy-MM-ddThh:mm:ss"));
      exp.setMetaValue("mzml_start_time_stamp", date);
    }
    copy_values(metadata.at("sample"), exp.getSample(), "Thermo ");
    exp.getSample().setName(text(metadata.at("sample"), "sample name"));
    exp.getSample().setNumber(text(metadata.at("sample"), "sample number"));
    exp.getSample().setComment(text(metadata.at("sample"), "sample comment"));
    exp.getSample().setMetaValue("Thermo user sample fields", metadata.at("user_sample_fields").dump());
    if (options_.instrument_methods) { exp.setMetaValue("Thermo instrument methods", metadata.at("instrument_methods").dump()); }
    if (metadata.at("run").is_object()) { copy_values(metadata.at("run"), exp, "Thermo run "); }
    exp.setMetaValue("Thermo RawFileReader version", text(metadata, "reader_version"));

    Instrument instrument;
    std::string model;
    if (metadata.at("instrument").is_object())
    {
      const auto& info = metadata.at("instrument");
      model = text(info, "model");
      instrument.setName(openms::thermo_bridge::cv_instrument_model(model).name);
      instrument.setModel(model);
      instrument.setMetaValue("Thermo instrument name", text(info, "name"));
      instrument.setMetaValue("instrument serial number", text(info, "serial_number"));
      instrument.setMetaValue("Thermo hardware version", text(info, "hardware_version"));
      instrument.setMetaValue("Thermo instrument units", text(info, "units"));
      instrument.getSoftware().setName("Thermo acquisition software");
      instrument.getSoftware().setVersion(text(info, "software_version"));
    }
    exp.setInstrument(instrument);
    auto processing = std::make_shared<DataProcessing>();
    processing->getSoftware().setName("OpenMS Thermo RAW reader");
    processing->getSoftware().setVersion("0.3.0");
    processing->setMetaValue("Thermo RawFileReader version", text(metadata, "reader_version"));
    processing->getProcessingActions().insert(DataProcessing::ProcessingAction::FORMAT_CONVERSION);
    if (options_.centroid) { processing->getProcessingActions().insert(DataProcessing::ProcessingAction::PEAK_PICKING); }

    std::map<std::string, std::string> configurations;
    Metadata precursor_parser;
    std::map<int, Size> spectrum_by_scan;
    auto activation = [&](Precursor& precursor, const Json& reaction, bool supplemental) {
      const std::string type = text(reaction, "activation");
      if (reaction.value("collision_energy_valid", false) && reaction.at("collision_energy").is_number())
      {
        precursor.setMetaValue(supplemental ? "supplemental collision energy" : "collision energy",
                               unit_value(numeric(reaction, "collision_energy"), 266));
      }
      if (supplemental)
      {
        if (type == "HigherEnergyCollisionalDissociation") { precursor.setMetaValue("supplemental beam-type collision-induced dissociation", ""); }
        else if (type == "CollisionInducedDissociation") { precursor.setMetaValue("supplemental collision-induced dissociation", ""); }
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
    };
    if (raw.has_ms_data())
    {
      const int first = raw.first_scan_number(), last = raw.last_scan_number();
      exp.reserve(raw.scan_count());
      startProgress(first, last, "Loading Thermo RAW file");
      for (int scan = first; scan <= last; ++scan)
      {
        setProgress(scan);
        const Json meta = Json::parse(raw.scan_metadata_json(scan));
        if (meta.at("schema_version") != 1) { throw std::runtime_error("Unsupported Thermo scan metadata schema"); }
        MSSpectrum spectrum;
        const int level = meta.at("ms_level");
        spectrum.setMSLevel(level > 0 ? level : 0);
        const std::string order = text(meta, "ms_order_name");
        auto mode = level == 1 ? InstrumentSettings::ScanMode::MS1SPECTRUM : InstrumentSettings::ScanMode::MSNSPECTRUM;
        if (order == "Par") { mode = InstrumentSettings::ScanMode::PRECURSOR; }
        else if (order == "Nl") { mode = InstrumentSettings::ScanMode::CNL; }
        else if (order == "Ng") { mode = InstrumentSettings::ScanMode::CNG; }
        spectrum.getInstrumentSettings().setScanMode(mode);
        spectrum.setRT(meta.at("retention_time").get<double>() * 60.0);
        spectrum.setNativeID(text(meta, "native_id"));
        const bool centroid = options_.centroid || meta.at("centroid").get<bool>();
        spectrum.setType(centroid ? SpectrumSettings::SpectrumType::CENTROID : SpectrumSettings::SpectrumType::PROFILE);
        const int polarity = meta.at("polarity");
        if (polarity == 1) { spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE); }
        else if (polarity == 0) { spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE); }
        const std::string filter = text(meta, "filter");
        Acquisition acquisition;
        acquisition.setMetaValue("filter string", filter);
        // Keep the historical spectrum-level accessor as well as the mzML
        // scan-level term.
        spectrum.setMetaValue("filter string", filter);
        if (options_.preserve_trailers) { acquisition.setMetaValue("Thermo trailer extra", meta.at("trailer").dump()); }
        auto injection = Metadata::number(Metadata::trailer(meta, "Ion Injection Time (ms):"));
        if (injection.is_number()) { acquisition.setMetaValue("ion injection time", unit_value(injection.get<double>(), 28)); }
        auto mono = Metadata::number(Metadata::trailer(meta, "Monoisotopic M/Z:"));
        if (mono.is_number() && mono.get<double>() > 0) { acquisition.setMetaValue("[Thermo Trailer Extra]Monoisotopic M/Z:", mono.get<double>()); }
        const auto voltage_on = Metadata::number(Metadata::trailer(meta, "FAIMS Voltage On:"));
        const auto voltage = Metadata::number(Metadata::trailer(meta, "FAIMS CV:"));
        std::string enabled = Metadata::trailer(meta, "FAIMS Voltage On:");
        std::transform(enabled.begin(), enabled.end(), enabled.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (((voltage_on.is_number() && voltage_on.get<double>() != 0) || (enabled == "true" || enabled == "on" || enabled == "yes"))
            && voltage.is_number())
        {
          spectrum.setDriftTime(voltage.get<double>());
          spectrum.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
        }
        const std::string analyzer = text(meta, "analyzer"), ionization = text(meta, "ionization");
        const std::string configuration_key = analyzer + ":" + ionization;
        if (! configurations.count(configuration_key))
        {
          const std::string id = configurations.empty() ? "ic_0" : "thermo_ic_" + std::to_string(configurations.size());
          configurations[configuration_key] = id;
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
          if (id != "ic_0") { exp.getInstrumentConfigurations()[id] = config; }
          if (configurations.size() == 1) { exp.setInstrument(config); }
        }
        if (configurations.at(configuration_key) != "ic_0")
        {
          acquisition.setMetaValue("instrument_configuration_ref", configurations.at(configuration_key));
        }
        spectrum.getAcquisitionInfo().push_back(acquisition);
        ScanWindow window;
        window.begin = numeric(meta, "low_mass");
        window.end = numeric(meta, "high_mass");
        if (meta.at("low_mass").is_number() && meta.at("high_mass").is_number() && window.end >= window.begin)
        {
          spectrum.getInstrumentSettings().getScanWindows().push_back(window);
        }
        if (meta.at("tic").is_number()) { spectrum.setMetaValue("total ion current", numeric(meta, "tic")); }
        if (meta.at("base_peak_mass").is_number())
        {
          spectrum.setMetaValue("base peak m/z", ms_unit_value(numeric(meta, "base_peak_mass"), 1000040));
        }
        if (meta.at("base_peak_intensity").is_number())
        {
          spectrum.setMetaValue("base peak intensity", ms_unit_value(numeric(meta, "base_peak_intensity"), 1000131));
        }
        const bool neutral = order == "Nl" || order == "Ng";
        if (neutral && ! meta.at("reactions").empty())
        {
          const auto& reaction = meta.at("reactions").front();
          spectrum.setMetaValue(order == "Nl" ? "neutral loss" : "neutral gain", numeric(reaction, "precursor_mass"));
          Product product;
          const double width = numeric(reaction, "isolation_width", -1);
          const double offset = numeric(reaction, "isolation_offset");
          if (width >= 0 && width / 2 >= std::abs(offset))
          {
            product.setIsolationWindowLowerOffset(width / 2 - offset);
            product.setIsolationWindowUpperOffset(width / 2 + offset);
          }
          spectrum.getProducts().push_back(product);
        }
        for (const auto& descriptor : neutral ? Json::array() : precursor_parser.precursors(meta))
        {
          Precursor precursor;
          precursor.setMZ(numeric(descriptor, "selected_mz")); // OpenMS default: selected ion
          precursor.setMetaValue("isolation window target m/z", numeric(descriptor, "target_mz"));
          precursor.setMetaValue("selected ion m/z", numeric(descriptor, "selected_mz"));
          if (descriptor.at("charge").is_number()) { precursor.setCharge(static_cast<int>(numeric(descriptor, "charge"))); }
          if (descriptor.at("width").is_number())
          {
            const double lower = numeric(descriptor, "lower_offset"), upper = numeric(descriptor, "upper_offset");
            if (lower >= 0) { precursor.setIsolationWindowLowerOffset(lower); }
            if (upper >= 0) { precursor.setIsolationWindowUpperOffset(upper); }
            if (lower <= 0) { precursor.setMetaValue("isolation window lower offset", lower); }
            if (upper <= 0) { precursor.setMetaValue("isolation window upper offset", upper); }
          }
          if (! text(descriptor, "spectrum_ref").empty()) { precursor.setMetaValue("spectrum_ref", text(descriptor, "spectrum_ref")); }
          activation(precursor, descriptor.at("activation"), false);
          if (descriptor.at("supplemental").is_object()) { activation(precursor, descriptor.at("supplemental"), true); }
          auto parent = spectrum_by_scan.find(descriptor.at("parent_scan").get<int>());
          if (descriptor.at("estimate_intensity").get<bool>() && parent != spectrum_by_scan.end() && precursor.getMZ() > 0)
          {
            // TRFP's precursor-intensity estimate sums the target +/- 1.5 m/z.
            double intensity = 0;
            const auto& parent_spectrum = exp[parent->second];
            const double target = numeric(descriptor, "target_mz");
            const double half_width = numeric(descriptor, "width") > 0 ? 1.5 : 0.0;
            for (auto peak = parent_spectrum.MZBegin(target - half_width); peak != parent_spectrum.end() && peak->getMZ() < target + half_width;
                 ++peak)
            {
              intensity += peak->getIntensity();
            }
            precursor.setIntensity(intensity);
            if (intensity == 0) { precursor.setMetaValue("peak intensity", 0.0); }
            precursor.setMetaValue("peak intensity unit accession", "MS:1000131");
          }
          spectrum.getPrecursors().push_back(precursor);
        }
        auto data = raw.spectrum_data(scan, options_.centroid);
        spectrum.resize(data.mz.size());
        for (Size i = 0; i < data.mz.size(); ++i)
        {
          spectrum[i].setMZ(data.mz[i]);
          spectrum[i].setIntensity(data.intensities[i]);
        }
        if (! data.mz.empty())
        {
          if (centroid)
          {
            const Size basepeak = std::distance(data.intensities.begin(), std::max_element(data.intensities.begin(), data.intensities.end()));
            spectrum.setMetaValue("base peak m/z", ms_unit_value(data.mz[basepeak], 1000040));
            spectrum.setMetaValue("base peak intensity", ms_unit_value(data.intensities[basepeak], 1000131));
          }
          const auto [lo, hi] = std::minmax_element(data.mz.begin(), data.mz.end());
          spectrum.setMetaValue("lowest observed m/z", ms_unit_value(*lo, 1000040));
          spectrum.setMetaValue("highest observed m/z", ms_unit_value(*hi, 1000040));
        }
        if (options_.charge_data && options_.centroid)
        {
          auto charges = raw.spectrum_auxiliary_array(scan, 0);
          if (! charges.empty() && charges.size() == spectrum.size())
          {
            MSSpectrum::IntegerDataArray array;
            array.setName("charge array");
            for (double charge : charges)
            {
              array.push_back(static_cast<Int64>(charge));
            }
            spectrum.getIntegerDataArrays().push_back(std::move(array));
          }
        }
        // The independently sampled noise grid must not be sorted/filtered with
        // peaks. Store double lists; MzMLHandler serializes these as 64-bit
        // binary arrays.
        if (options_.noise_data)
        {
          const char* names[] = {"", "sampled noise m/z array", "sampled noise intensity array", "sampled noise baseline array"};
          for (int kind = 1; kind <= 3; ++kind)
          {
            auto values = raw.spectrum_auxiliary_array(scan, kind);
            if (! values.empty()) { spectrum.setMetaValue(names[kind], values); }
          }
        }
        spectrum.getDataProcessing().push_back(processing);
        if (! spectrum.isSorted()) { spectrum.sortByPosition(); }
        spectrum_by_scan[scan] = exp.size();
        exp.addSpectrum(std::move(spectrum));
      }
      endProgress();
      const auto data = raw.chromatogram_data();
      MSChromatogram tic;
      tic.setNativeID("TIC");
      tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
      for (Size i = 0; i < data.times.size(); ++i)
      {
        tic.push_back(ChromatogramPeak(data.times[i] * 60.0, data.intensities[i]));
      }
      tic.getDataProcessing().push_back(processing);
      if (! tic.empty()) { exp.addChromatogram(std::move(tic)); }
    }
    if (options_.all_detectors)
    {
      const Json detectors = Json::parse(raw.detector_chromatograms_json());
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
        for (Size i = 0; i < times.size(); ++i)
        {
          chromatogram.push_back(ChromatogramPeak(times[i] * 60.0, values[i]));
        }
        chromatogram.getDataProcessing().push_back(processing);
        exp.addChromatogram(std::move(chromatogram));
      }
      for (const auto& controller : metadata.at("controllers"))
      {
        if (text(controller, "name") != "Pda") { continue; }
        raw.select_instrument(controller.at("type").get<int>(), controller.at("number").get<int>());
        for (int scan = raw.first_scan_number(); scan <= raw.last_scan_number(); ++scan)
        {
          const Json meta = Json::parse(raw.scan_metadata_json(scan));
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
          if (meta.at("low_wavelength").is_number() && meta.at("high_wavelength").is_number())
          {
            spectrum.getInstrumentSettings().getScanWindows().push_back(window);
          }
          const auto data = raw.spectrum_data(scan, false);
          for (Size i = 0; i < data.mz.size(); ++i)
          {
            Peak1D peak;
            peak.setMZ(data.mz[i]);
            peak.setIntensity(data.intensities[i]);
            spectrum.push_back(peak);
          }
          if (! data.mz.empty())
          {
            const auto [lo, hi] = std::minmax_element(data.mz.begin(), data.mz.end());
            spectrum.setMetaValue("lowest observed wavelength", unit_value(*lo, 18));
            spectrum.setMetaValue("highest observed wavelength", unit_value(*hi, 18));
          }
          spectrum.getDataProcessing().push_back(processing);
          exp.addSpectrum(std::move(spectrum));
        }
      }
    }
    exp.sortSpectra(true);
    exp.updateRanges();
  }
  catch (const std::exception& error)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, std::string("Thermo RAW reader: ") + error.what());
  }
}

} // namespace OpenMS
#endif // WITH_THERMO_RAW
