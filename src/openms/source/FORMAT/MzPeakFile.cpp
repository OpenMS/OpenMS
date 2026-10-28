// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzPeakFile.h>

#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/IONMOBILITY/IMDataArrayUtils.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>

#ifdef WITH_MZPEAK
#include <mzpeak/chromatograms.h>
#include <mzpeak/exception.h>
#include <mzpeak/index.h>
#include <mzpeak/open.h>
#include <mzpeak/spectra.h>
#include <mzpeak/spectrum.h>
#include <mzpeak/writer.h>

#include <boost/json.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif
#endif

namespace OpenMS
{
  MzPeakFile::MzPeakFile() = default;

  MzPeakFile::~MzPeakFile() = default;

  PeakFileOptions& MzPeakFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& MzPeakFile::getOptions() const
  {
    return options_;
  }

  void MzPeakFile::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  bool MzPeakFile::isMzPeak(const std::string& path)
  {
    std::error_code ec;
    const std::filesystem::path p(path);
    if (std::filesystem::is_directory(p, ec))
    {
      return std::filesystem::is_regular_file(p / "mzpeak_index.json", ec);
    }
    if (!std::filesystem::is_regular_file(p, ec)) { return false; }
    std::ifstream in(path, std::ios::binary);
    char magic[4] = {0, 0, 0, 0};
    in.read(magic, 4);
    // a ZIP archive: local file header, or an empty archive's end-of-central-directory record
    return in.gcount() == 4 && magic[0] == 'P' && magic[1] == 'K' && ((magic[2] == 3 && magic[3] == 4) || (magic[2] == 5 && magic[3] == 6));
  }

  bool MzPeakFile::isSupported()
  {
#ifdef WITH_MZPEAK
    return true;
#else
    return false;
#endif
  }

#ifndef WITH_MZPEAK

  void MzPeakFile::load(const std::string& /*filename*/, PeakMap& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void MzPeakFile::transform(const std::string& /*filename*/, Interfaces::IMSDataConsumer* /*consumer*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void MzPeakFile::store(const std::string& /*filename*/, const PeakMap& /*exp*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

#else

  namespace
  {
    // ---------------------------------------------------------------------------------------------------------------
    // Reading
    // ---------------------------------------------------------------------------------------------------------------

    /// The spectra decoded in parallel per call of the consumer (transform()): bounds the memory of a streamed read
    constexpr Size BATCH_SIZE = 2048;

    /// How often the progress is updated (spectra)
    constexpr Size PROGRESS_STEP = 1000;

    /// Dissociation method CV accessions -> OpenMS activation methods (as MzMLHandler maps them)
    bool activationMethod(const std::string& accession, Precursor::ActivationMethod& method)
    {
      static const std::map<std::string, Precursor::ActivationMethod> methods = {
        {"MS:1000133", Precursor::ActivationMethod::CID},      {"MS:1000134", Precursor::ActivationMethod::PD},
        {"MS:1000135", Precursor::ActivationMethod::PSD},      {"MS:1000136", Precursor::ActivationMethod::SID},
        {"MS:1000242", Precursor::ActivationMethod::BIRD},     {"MS:1000250", Precursor::ActivationMethod::ECD},
        {"MS:1000262", Precursor::ActivationMethod::IMD},      {"MS:1000282", Precursor::ActivationMethod::SORI},
        {"MS:1000422", Precursor::ActivationMethod::HCD},      {"MS:1002472", Precursor::ActivationMethod::TRAP},
        {"MS:1002481", Precursor::ActivationMethod::HCID},     {"MS:1000433", Precursor::ActivationMethod::LCID},
        {"MS:1000435", Precursor::ActivationMethod::PHD},      {"MS:1000598", Precursor::ActivationMethod::ETD},
        {"MS:1002679", Precursor::ActivationMethod::ETciD},    {"MS:1003182", Precursor::ActivationMethod::ETciD},
        {"MS:1002678", Precursor::ActivationMethod::EThcD},    {"MS:1002631", Precursor::ActivationMethod::EThcD},
        {"MS:1000599", Precursor::ActivationMethod::PQD},      {"MS:1001880", Precursor::ActivationMethod::INSOURCE},
      };
      const auto it = methods.find(accession);
      if (it == methods.end()) { return false; }
      method = it->second;
      return true;
    }

    /// The OpenMS unit of an ion mobility CV accession (ion_mobility_type); NONE if unknown
    DriftTimeUnit driftTimeUnit(const std::optional<std::string>& accession)
    {
      if (!accession) { return DriftTimeUnit::NONE; }
      if (*accession == "MS:1002815") { return DriftTimeUnit::VSSC; }                       // inverse reduced ion mobility
      if (*accession == "MS:1002476") { return DriftTimeUnit::MILLISECOND; }                // ion mobility drift time
      if (*accession == "MS:1001581") { return DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE; } // FAIMS compensation voltage
      if (*accession == "MS:1002954") { return DriftTimeUnit::CCS; }                        // collisional cross sectional area
      return DriftTimeUnit::NONE;
    }

    /// The CV accession of an OpenMS drift time unit (for writing); empty if none
    std::optional<std::string> ionMobilityType(DriftTimeUnit unit)
    {
      switch (unit)
      {
        case DriftTimeUnit::VSSC: return "MS:1002815";
        case DriftTimeUnit::MILLISECOND: return "MS:1002476";
        case DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE: return "MS:1001581";
        case DriftTimeUnit::CCS: return "MS:1002954";
        default: return std::nullopt;
      }
    }

    /// A CV parameter's value as a DataValue: a number if it parses as one, else the string
    DataValue cvValue(const MzPeak::CvParam& p)
    {
      if (!p.value) { return DataValue(); }
      const std::string& v = *p.value;
      if (v.empty()) { return DataValue(v); }
      char* end = nullptr;
      const double d = std::strtod(v.c_str(), &end);
      if (end != nullptr && *end == '\0' && end != v.c_str())
      {
        if (v.find_first_of(".eEnN") == std::string::npos && d >= -2147483648.0 && d <= 2147483647.0) { return DataValue(static_cast<Int>(d)); }
        return DataValue(d);
      }
      return DataValue(v);
    }

    /// Store a CV parameter as a meta value named after the term (as MzMLHandler does for terms without a member)
    void setCvMeta(MetaInfoInterface& target, const MzPeak::CvParam& p)
    {
      if (!p.name || p.name->empty() || !p.value) { return; }
      target.setMetaValue(*p.name, cvValue(p));
    }

    /// Whether the spectrum passes the scan-level options (MS level, RT range, precursor m/z range), from its metadata
    /// alone (no peak decode)
    bool keepByMetadata(const MzPeak::SpectrumMetadata& md, const PeakFileOptions& options)
    {
      if (options.hasMSLevels())
      {
        if (!md.ms_level || !options.containsMSLevel(*md.ms_level)) { return false; }
      }
      if (options.hasRTRange())
      {
        if (!md.retention_time || !options.getRTRange().encloses(DPosition<1>(*md.retention_time))) { return false; }
      }
      if (options.hasPrecursorMZRange() && md.ms_level.value_or(1) > 1 && !md.precursors.empty())
      {
        const MzPeak::PrecursorInfo& p = md.precursors.front();
        std::optional<double> mz;
        if (!p.selected_ions.empty() && p.selected_ions.front().selected_ion_mz && options.getPrecursorMZSelectedIon())
        {
          mz = p.selected_ions.front().selected_ion_mz;
        }
        else if (p.isolation_window.target_mz)
        {
          mz = *p.isolation_window.target_mz;
        }
        else if (!p.selected_ions.empty() && p.selected_ions.front().selected_ion_mz)
        {
          mz = p.selected_ions.front().selected_ion_mz;
        }
        if (mz && !options.getPrecursorMZRange().encloses(DPosition<1>(*mz))) { return false; }
      }
      return true;
    }

    /// The precursor of an mzPeak precursor record (isolation window, first selected ion, activation), as MzMLHandler
    /// builds it from mzML
    Precursor toPrecursor(const MzPeak::PrecursorInfo& in, const PeakFileOptions& options)
    {
      Precursor out;
      const MzPeak::IsolationWindow& iw = in.isolation_window;
      bool have_target = false;
      if (iw.target_mz)
      {
        out.setMZ(*iw.target_mz);
        have_target = true;
      }
      if (iw.lower_offset && *iw.lower_offset >= 0) { out.setIsolationWindowLowerOffset(*iw.lower_offset); }
      if (iw.upper_offset && *iw.upper_offset >= 0) { out.setIsolationWindowUpperOffset(*iw.upper_offset); }
      for (const MzPeak::CvParam& p : iw.parameters) { setCvMeta(out, p); }
      if (!in.selected_ions.empty())
      {
        const MzPeak::SelectedIonInfo& ion = in.selected_ions.front();
        if (ion.selected_ion_mz)
        {
          const double mz = *ion.selected_ion_mz;
          // The isolation target is a float in mzPeak: it is the same value as the selected ion m/z if they agree at float
          // precision, and then (as for mzML, where the two are usually the same number) neither meta value is written.
          const bool same = have_target && static_cast<float>(mz) == static_cast<float>(out.getMZ());
          if (!have_target || options.getPrecursorMZSelectedIon())
          {
            if (have_target && !same) { out.setMetaValue("isolation window target m/z", out.getMZ()); }
            out.setMZ(mz); // (the double, also when the float target agrees)
          }
          else if (!same)
          {
            out.setMetaValue("selected ion m/z", mz);
          }
        }
        if (ion.charge_state) { out.setCharge(*ion.charge_state); }
        if (ion.intensity) { out.setIntensity(*ion.intensity); }
        if (ion.ion_mobility_value)
        {
          out.setDriftTime(*ion.ion_mobility_value);
          out.setDriftTimeUnit(driftTimeUnit(ion.ion_mobility_type));
        }
        for (const MzPeak::CvParam& p : ion.parameters) { setCvMeta(out, p); }
        if (in.selected_ions.size() > 1) { out.setMetaValue("selected ions", static_cast<Int>(in.selected_ions.size())); }
      }
      for (const MzPeak::CvParam& p : in.activation_parameters)
      {
        if (!p.accession) { continue; }
        Precursor::ActivationMethod method;
        if (activationMethod(*p.accession, method))
        {
          out.getActivationMethods().insert(method);
        }
        else if (*p.accession == "MS:1000045" && p.value) // collision energy (eV)
        {
          out.setActivationEnergy(std::strtod(p.value->c_str(), nullptr));
        }
        else if (*p.accession != "MS:1000044") // (dissociation method: the parent term carries no information)
        {
          setCvMeta(out, p);
        }
      }
      if (!in.precursor_id.empty()) { out.setMetaValue("spectrum_ref", in.precursor_id); }
      return out;
    }

    /// Keep only the peaks of @p spectrum within the m/z, intensity and ion mobility ranges of @p options (the per-peak
    /// data arrays stay aligned). Returns false if the spectrum itself is excluded (a single ion mobility outside the range).
    bool filterPeaks(const PeakFileOptions& options, MSSpectrum& spectrum)
    {
      if (options.hasIMRange() && !spectrum.containsIMData() && spectrum.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET
          && !options.getIMRange().encloses(DPosition<1>(spectrum.getDriftTime())))
      {
        return false;
      }
      const bool im = options.hasIMRange() && spectrum.containsIMData();
      if (!options.hasMZRange() && !options.hasIntensityRange() && !im) { return true; }
      const DataArrays::FloatDataArray* im_array = im ? &spectrum.getFloatDataArrays()[spectrum.getIMData().first] : nullptr;
      std::vector<Size> kept;
      kept.reserve(spectrum.size());
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        if (options.hasMZRange() && !options.getMZRange().encloses(DPosition<1>(spectrum[i].getMZ()))) { continue; }
        if (options.hasIntensityRange() && !options.getIntensityRange().encloses(DPosition<1>(spectrum[i].getIntensity()))) { continue; }
        if (im_array != nullptr && i < im_array->size() && !options.getIMRange().encloses(DPosition<1>((*im_array)[i]))) { continue; }
        kept.push_back(i);
      }
      if (kept.size() != spectrum.size()) { spectrum.select(kept); }
      return true;
    }

    /// The MSSpectrum of an mzPeak spectrum: metadata, and the peaks unless @p metadata_only
    void toSpectrum(const MzPeak::Spectrum& in, MSSpectrum& out, const PeakFileOptions& options, bool metadata_only)
    {
      const MzPeak::SpectrumMetadata& md = in.metadata();
      out.setNativeID(md.id);
      if (md.ms_level) { out.setMSLevel(static_cast<UInt>(*md.ms_level)); }
      if (md.retention_time) { out.setRT(*md.retention_time); } // seconds (the library converts from the minutes stored)
      if (md.polarity)
      {
        if (*md.polarity > 0) { out.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE); }
        else if (*md.polarity < 0) { out.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE); }
      }
      if (md.representation == "MS:1000127") { out.setType(SpectrumSettings::SpectrumType::CENTROID); }
      else if (md.representation == "MS:1000128") { out.setType(SpectrumSettings::SpectrumType::PROFILE); }
      if (md.base_peak_mz) { out.setMetaValue("base peak m/z", *md.base_peak_mz); }
      if (md.base_peak_intensity) { out.setMetaValue("base peak intensity", static_cast<double>(*md.base_peak_intensity)); }
      if (md.total_ion_current) { out.setMetaValue("total ion current", static_cast<double>(*md.total_ion_current)); }
      if (md.lowest_observed_mz) { out.setMetaValue("lowest observed m/z", *md.lowest_observed_mz); }
      if (md.highest_observed_mz) { out.setMetaValue("highest observed m/z", *md.highest_observed_mz); }
      if (!md.spectrum_type.empty()) { out.setMetaValue("spectrum type", md.spectrum_type); }
      for (const MzPeak::CvParam& p : md.parameters)
      {
        if (p.accession && *p.accession == "MS:1001581" && p.value) // FAIMS compensation voltage (as a spectrum term)
        {
          out.setDriftTime(std::strtod(p.value->c_str(), nullptr));
          out.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
          continue;
        }
        setCvMeta(out, p);
      }
      for (const MzPeak::CvParam& p : md.scan_parameters)
      {
        if (p.accession && *p.accession == "MS:1000016") { continue; } // scan start time: the RT
        if (p.accession && (*p.accession == "MS:1002815" || *p.accession == "MS:1002476" || *p.accession == "MS:1001581" || *p.accession == "MS:1002954") && p.value)
        {
          out.setDriftTime(std::strtod(p.value->c_str(), nullptr));
          out.setDriftTimeUnit(driftTimeUnit(p.accession));
          continue;
        }
        setCvMeta(out, p);
      }
      if (md.ion_mobility) // the scan's ion mobility (ion_mobility_value / ion_mobility_type of the scans facet)
      {
        out.setDriftTime(*md.ion_mobility);
        out.setDriftTimeUnit(driftTimeUnit(md.ion_mobility_type));
      }
      for (const MzPeak::ScanWindow& w : md.scan_windows)
      {
        ScanWindow sw;
        if (w.lower_limit) { sw.begin = *w.lower_limit; }
        if (w.upper_limit) { sw.end = *w.upper_limit; }
        for (const MzPeak::CvParam& p : w.parameters) { setCvMeta(sw, p); }
        out.getInstrumentSettings().getScanWindows().push_back(sw);
      }
      for (const MzPeak::PrecursorInfo& p : md.precursors) { out.getPrecursors().push_back(toPrecursor(p, options)); }
      if (metadata_only) { return; }

      const std::vector<double>& mz = in.mz();
      const std::vector<float>& intensity = in.intensity();
      const Size n = std::min(mz.size(), intensity.size());
      if (mz.size() != intensity.size())
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, md.id,
                                    "mzPeak spectrum with " + std::to_string(mz.size()) + " m/z and " + std::to_string(intensity.size()) + " intensity values");
      }
      out.resize(n);
      for (Size i = 0; i < n; ++i)
      {
        out[i].setMZ(mz[i]);
        out[i].setIntensity(intensity[i]);
      }
      const std::vector<double>& mobility = in.ion_mobility_array();
      if (!mobility.empty() && mobility.size() == n)
      {
        DataArrays::FloatDataArray im_array;
        im_array.resize(n);
        for (Size i = 0; i < n; ++i) { im_array[i] = static_cast<float>(mobility[i]); }
        DriftTimeUnit unit = driftTimeUnit(md.ion_mobility_type);
        if (unit == DriftTimeUnit::NONE && !md.precursors.empty() && !md.precursors.front().selected_ions.empty())
        {
          unit = driftTimeUnit(md.precursors.front().selected_ions.front().ion_mobility_type);
        }
        if (unit == DriftTimeUnit::NONE) { unit = DriftTimeUnit::VSSC; } // (timsTOF data: 1/K0 is what the converters write)
        IMDataArrayUtils::setIMUnit(im_array, unit);
        out.getFloatDataArrays().push_back(std::move(im_array));
        if (md.ion_mobility) // a spectrum with an ion mobility per peak holds no single drift time
        {
          out.setDriftTime(IMTypes::DRIFTTIME_NOT_SET);
          out.setDriftTimeUnit(DriftTimeUnit::NONE);
        }
      }
    }

    /// The chromatogram type of the mzPeak chromatogram type accession (as MzMLHandler maps them); false if unknown
    bool chromatogramType(const std::string& accession, ChromatogramSettings::ChromatogramType& type)
    {
      using CT = ChromatogramSettings::ChromatogramType;
      static const std::map<std::string, CT> types = {
        {"MS:1000810", CT::MASS_CHROMATOGRAM},
        {"MS:1000235", CT::TOTAL_ION_CURRENT_CHROMATOGRAM},
        {"MS:1000627", CT::SELECTED_ION_CURRENT_CHROMATOGRAM},
        {"MS:1000628", CT::BASEPEAK_CHROMATOGRAM},
        {"MS:1001472", CT::SELECTED_ION_MONITORING_CHROMATOGRAM},
        {"MS:1001473", CT::SELECTED_REACTION_MONITORING_CHROMATOGRAM},
        {"MS:1001474", CT::SELECTED_REACTION_MONITORING_CHROMATOGRAM},
        {"MS:1000811", CT::ELECTROMAGNETIC_RADIATION_CHROMATOGRAM},
        {"MS:1000812", CT::ABSORPTION_CHROMATOGRAM},
        {"MS:1000813", CT::EMISSION_CHROMATOGRAM},
      };
      const auto it = types.find(accession);
      if (it == types.end()) { return false; }
      type = it->second;
      return true;
    }

    /// The MSChromatogram of an mzPeak chromatogram
    MSChromatogram toChromatogram(const MzPeak::Chromatogram& in, const PeakFileOptions& options)
    {
      MSChromatogram out;
      const MzPeak::ChromatogramMetadata& md = in.metadata();
      out.setNativeID(md.id);
      ChromatogramSettings::ChromatogramType type;
      if (chromatogramType(md.chromatogram_type, type)) { out.setChromatogramType(type); }
      else if (!md.chromatogram_type.empty()) { out.setMetaValue("chromatogram type accession", md.chromatogram_type); }
      for (const MzPeak::CvParam& p : md.parameters) { setCvMeta(out, p); }
      if (!md.precursors.empty()) { out.setPrecursor(toPrecursor(md.precursors.front(), options)); }
      const std::vector<double>& time = in.time(); // seconds
      const std::vector<float>& intensity = in.intensity();
      const Size n = std::min(time.size(), intensity.size());
      out.resize(n);
      for (Size i = 0; i < n; ++i)
      {
        out[i].setRT(time[i]);
        out[i].setIntensity(intensity[i]);
      }
      return out;
    }

    /// The experimental settings of the run-level metadata of the archive
    void toSettings(const MzPeak::Index& index, const std::string& filename, ExperimentalSettings& settings)
    {
      const MzPeak::RunMetadata& md = index.metadata();
      if (md.run())
      {
        if (md.run()->id) { settings.setIdentifier(*md.run()->id); }
        if (md.run()->start_time)
        {
          try
          {
            std::string t = *md.run()->start_time;
            if (t.size() > 19) { t.resize(19); } // "yyyy-MM-ddThh:mm:ss" (no fraction, no zone)
            settings.setDateTime(DateTime::fromString(t));
          }
          catch (...)
          {
            settings.setMetaValue("run start time", *md.run()->start_time);
          }
        }
      }
      if (md.file_description())
      {
        for (const MzPeak::SourceFile& sf : md.file_description()->source_files)
        {
          SourceFile out;
          if (sf.name) { out.setNameOfFile(*sf.name); }
          if (sf.location) { out.setPathToFile(*sf.location); }
          for (const MzPeak::CvParam& p : sf.parameters)
          {
            if (!p.accession || !p.name) { continue; }
            if (p.accession->rfind("MS:", 0) == 0 && p.name->find("file format") != std::string::npos) { out.setFileType(*p.name); }
            else if (*p.accession == "MS:1000569" && p.value) { out.setChecksum(*p.value, SourceFile::ChecksumType::SHA1); }
            else if (*p.accession == "MS:1000568" && p.value) { out.setChecksum(*p.value, SourceFile::ChecksumType::MD5); }
            else { setCvMeta(out, p); }
          }
          settings.getSourceFiles().push_back(out);
        }
        for (const MzPeak::CvParam& p : md.file_description()->contents) { setCvMeta(settings, p); }
      }
      if (!md.instrument_configurations().empty())
      {
        const MzPeak::InstrumentConfiguration& ic = md.instrument_configurations().front();
        for (const MzPeak::CvParam& p : ic.parameters)
        {
          if (!p.accession || !p.name) { continue; }
          if (*p.accession == "MS:1000529" && p.value) { settings.getInstrument().setMetaValue("instrument serial number", *p.value); }
          else if (!p.value || p.value->empty()) { settings.getInstrument().setName(*p.name); } // the instrument model term
          else { setCvMeta(settings.getInstrument(), p); }
        }
      }
      // (the software list has no place in ExperimentalSettings; data processing is per spectrum in OpenMS)
      (void)filename;
    }

    /// Opens the archive, with the library's errors as OpenMS exceptions
    MzPeak::Index openIndex(const std::string& filename)
    {
      if (!File::exists(filename)) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename); }
      try
      {
        return MzPeak::open(filename);
      }
      catch (const std::exception& e)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, std::string("cannot open mzPeak run: ") + e.what());
      }
    }

    /// The spectra of an archive, read with the OpenMP threads of the calling thread: one Spectra (the library's
    /// reader) per thread, the metadata shared through the archive
    class SpectraReader
    {
    public:
      SpectraReader(const MzPeak::Index& index, const PeakFileOptions& options) :
        index_(index), options_(options), spectra_(index.spectra())
      {
        // the spectra the scan-level options keep, from the metadata alone
        const Size n = spectra_.size();
        kept_.reserve(n);
        for (Size i = 0; i < n; ++i)
        {
          if (keepByMetadata(spectra_[i].metadata(), options_)) { kept_.push_back(i); }
        }
      }

      /// The number of spectra in the file
      Size total() const { return spectra_.size(); }

      /// The indices (in file order) of the spectra the options keep
      const std::vector<Size>& kept() const { return kept_; }

      /// The spectra kept[first..last) as MSSpectra, decoded in parallel: out[k] is the spectrum of kept[first + k];
      /// excluded[k] is set for a spectrum the peak-level options exclude (left empty).
      void decode(Size first, Size last, std::vector<MSSpectrum>& out, std::vector<char>& excluded, const std::string& filename)
      {
        out.clear();
        out.resize(last - first);
        excluded.assign(last - first, 0);
        std::exception_ptr error;
        const bool metadata_only = options_.getMetadataOnly() || !options_.getFillData();
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
          // one reader per thread (the library's concurrency model); the archive's metadata map is shared
          std::unique_ptr<MzPeak::Spectra> local;
#ifdef _OPENMP
#pragma omp critical(MzPeakFile_reader)
#endif
          {
            try
            {
              local.reset(new MzPeak::Spectra(index_.spectra())); // (Spectra is not movable: direct-initialized from the prvalue)
            }
            catch (...)
            {
              if (!error) { error = std::current_exception(); }
            }
          }
          if (local)
          {
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
            for (SignedSize k = 0; k < static_cast<SignedSize>(last - first); ++k)
            {
              try
              {
                if (error) { continue; }
                MSSpectrum& s = out[k];
                toSpectrum((*local)[kept_[first + k]], s, options_, metadata_only);
                if (!filterPeaks(options_, s))
                {
                  s = MSSpectrum();
                  excluded[k] = 1;
                }
              }
              catch (...)
              {
#ifdef _OPENMP
#pragma omp critical(MzPeakFile_error)
#endif
                {
                  if (!error) { error = std::current_exception(); }
                }
              }
            }
          }
        }
        if (error)
        {
          try
          {
            std::rethrow_exception(error);
          }
          catch (Exception::BaseException&)
          {
            throw;
          }
          catch (const std::exception& e)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, std::string("error reading mzPeak spectra: ") + e.what());
          }
        }
      }

    private:
      const MzPeak::Index& index_;
      const PeakFileOptions& options_;
      MzPeak::Spectra spectra_;
      std::vector<Size> kept_;
    };

    /// The chromatograms of an archive
    std::vector<MSChromatogram> readChromatograms(const MzPeak::Index& index, const PeakFileOptions& options, const std::string& filename)
    {
      std::vector<MSChromatogram> out;
      try
      {
        MzPeak::Chromatograms chromatograms = index.chromatograms();
        const Size n = chromatograms.size();
        out.reserve(n);
        for (Size i = 0; i < n; ++i) { out.push_back(toChromatogram(chromatograms[i], options)); }
      }
      catch (Exception::BaseException&)
      {
        throw;
      }
      catch (const std::exception& e)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, std::string("error reading mzPeak chromatograms: ") + e.what());
      }
      return out;
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Writing
    // ---------------------------------------------------------------------------------------------------------------

    MzPeak::SpectrumData toSpectrumData(const MSSpectrum& in, Size index)
    {
      MzPeak::SpectrumData out;
      const Size n = in.size();
      out.mz.resize(n);
      out.intensity.resize(n);
      for (Size i = 0; i < n; ++i)
      {
        out.mz[i] = in[i].getMZ();
        out.intensity[i] = in[i].getIntensity();
      }
      out.centroid = in.getType() != SpectrumSettings::SpectrumType::PROFILE; // (unknown: centroided, as most OpenMS inputs are)
      out.ms_level = static_cast<uint8_t>(std::min<UInt>(in.getMSLevel(), 255));
      out.retention_time = in.getRT();
      switch (in.getInstrumentSettings().getPolarity())
      {
        case IonSource::Polarity::POSITIVE: out.polarity = 1; break;
        case IonSource::Polarity::NEGATIVE: out.polarity = -1; break;
        default: break;
      }
      out.id = in.getNativeID().empty() ? "index=" + std::to_string(index) : in.getNativeID();
      for (const Precursor& p : in.getPrecursors())
      {
        MzPeak::PrecursorData pd;
        pd.isolation_target_mz = static_cast<float>(p.getMZ());
        pd.isolation_lower_offset = static_cast<float>(p.getIsolationWindowLowerOffset());
        pd.isolation_upper_offset = static_cast<float>(p.getIsolationWindowUpperOffset());
        MzPeak::SelectedIonData ion;
        ion.mz = p.getMZ();
        if (p.getCharge() != 0) { ion.charge = p.getCharge(); }
        if (p.getIntensity() > 0) { ion.intensity = p.getIntensity(); }
        pd.selected_ions.push_back(ion);
        out.precursors.push_back(std::move(pd));
      }
      if (in.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET)
      {
        out.ion_mobility = in.getDriftTime();
        out.ion_mobility_type = ionMobilityType(in.getDriftTimeUnit());
      }
      return out;
    }

    /// The run-level metadata of the archive: run id, source files, software
    MzPeak::RunMetadata toRunMetadata(const PeakMap& exp)
    {
      namespace json = boost::json;
      json::object root;
      json::object run;
      if (!exp.getIdentifier().empty()) { run["id"] = exp.getIdentifier(); }
      if (exp.getDateTime().isValid()) { run["start_time"] = exp.getDateTime().get(); }
      if (!run.empty()) { root["run"] = run; }
      json::array source_files;
      Size k = 0;
      for (const SourceFile& sf : exp.getSourceFiles())
      {
        json::object o;
        o["id"] = "sf" + std::to_string(k++);
        o["name"] = sf.getNameOfFile();
        o["location"] = sf.getPathToFile();
        source_files.push_back(o);
      }
      if (!source_files.empty())
      {
        json::object fd;
        fd["source_files"] = source_files;
        root["file_description"] = fd;
      }
      json::array software;
      {
        json::object o;
        o["id"] = "OpenMS";
        o["version"] = VersionInfo::getVersion();
        software.push_back(o);
      }
      root["software_list"] = software;
      return MzPeak::RunMetadata(root);
    }
  } // namespace

  void MzPeakFile::load(const std::string& filename, PeakMap& exp)
  {
    exp.reset();
    const MzPeak::Index index = openIndex(filename);
    toSettings(index, filename, exp);
    SpectraReader reader(index, options_);
    const std::vector<Size>& kept = reader.kept();
    startProgress(0, kept.size(), "Loading mzPeak file");
    std::vector<MSSpectrum> spectra;
    std::vector<char> excluded;
    reader.decode(0, kept.size(), spectra, excluded, filename);
    endProgress();
    // drop the spectra the peak-level options excluded
    Size n = 0;
    for (Size i = 0; i < spectra.size(); ++i)
    {
      if (excluded[i]) { continue; }
      if (n != i) { spectra[n] = std::move(spectra[i]); }
      ++n;
    }
    spectra.resize(n);
    exp.setSpectra(std::move(spectra));
    exp.setChromatograms(readChromatograms(index, options_, filename));
    exp.setLoadedFileType(FileTypes::MZPEAK); // (not by content: DocumentIdentifier would open the ZIP as a compressed single file)
    exp.setLoadedFilePath(filename);
    exp.updateRanges();
  }

  void MzPeakFile::transform(const std::string& filename, Interfaces::IMSDataConsumer* consumer)
  {
    if (consumer == nullptr) { throw Exception::NullPointer(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    const MzPeak::Index index = openIndex(filename);
    ExperimentalSettings settings;
    toSettings(index, filename, settings);
    SpectraReader reader(index, options_);
    const std::vector<Size>& kept = reader.kept();
    std::vector<MSChromatogram> chromatograms = readChromatograms(index, options_, filename);
    consumer->setExperimentalSettings(settings);
    consumer->setExpectedSize(kept.size(), chromatograms.size());
    startProgress(0, kept.size(), "Reading mzPeak file");
    std::vector<MSSpectrum> batch;
    std::vector<char> excluded;
    for (Size first = 0; first < kept.size(); first += BATCH_SIZE)
    {
      const Size last = std::min(first + BATCH_SIZE, kept.size());
      reader.decode(first, last, batch, excluded, filename);
      for (Size k = 0; k < batch.size(); ++k)
      {
        if (excluded[k]) { continue; } // by the peak-level options
        consumer->consumeSpectrum(batch[k]);
      }
      setProgress(last);
    }
    endProgress();
    for (MSChromatogram& c : chromatograms) { consumer->consumeChromatogram(c); }
  }

  void MzPeakFile::store(const std::string& filename, const PeakMap& exp) const
  {
    try
    {
      std::error_code ec;
      std::filesystem::remove(filename, ec); // libzip would add to an existing archive
      MzPeak::RunArchiveWriter writer(filename);
      writer.set_metadata(toRunMetadata(exp));
      startProgress(0, exp.size(), "Storing mzPeak file");
      for (Size i = 0; i < exp.size(); ++i)
      {
        writer.add(toSpectrumData(exp[i], i));
        if (i % PROGRESS_STEP == 0) { setProgress(i); }
      }
      writer.finish();
      endProgress();
    }
    catch (Exception::BaseException&)
    {
      throw;
    }
    catch (const std::exception& e)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename + " (" + e.what() + ")");
    }
    if (!exp.getChromatograms().empty())
    {
      OPENMS_LOG_WARN << "MzPeakFile: the " << exp.getChromatograms().size() << " chromatograms of the experiment were not written to '" << filename << "' (not supported by the writer yet)" << std::endl;
    }
  }

#endif // WITH_MZPEAK

} // namespace OpenMS
