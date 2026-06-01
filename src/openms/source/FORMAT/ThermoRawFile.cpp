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
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <filesystem>

#include <openms_thermo_bridge/thermo_bridge.hpp>
#include <openms_thermo_bridge/cv_mapping.hpp>

namespace OpenMS
{

  void ThermoRawFile::load(const String& path, MSExperiment& exp)
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

      // --- Instrument metadata ---
      {
        const std::string model = raw.instrument_model();
        if (!model.empty())
        {
          exp.getInstrument().setName(raw.instrument_name());
          exp.getInstrument().setModel(model);
          exp.getInstrument().setCustomizations(raw.instrument_serial_number());
        }
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

        // Polarity — Thermo PolarityType enum: Positive = 1, Negative = 2, Any = 0
        const int polarity = raw.polarity(scan);
        if (polarity == 1)
        {
          spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::POSITIVE);
        }
        else if (polarity == 2)
        {
          spectrum.getInstrumentSettings().setPolarity(IonSource::Polarity::NEGATIVE);
        }

        // Scan filter as comment
        const std::string filter = raw.scan_filter(scan);
        if (!filter.empty())
        {
          spectrum.setComment(filter);
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
          const double ce = raw.collision_energy(scan);
          if (ce > 0.0)
          {
            // Map activation type string to OpenMS activation method
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
            else if (!act_type.empty())
            {
              OPENMS_LOG_WARN << "ThermoRawFile: unknown activation type '" << act_type
                              << "' for scan " << scan << ", defaulting to CID\n";
              prec.getActivationMethods().insert(Precursor::ActivationMethod::CID);
            }
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

        // Populate peaks
        spectrum.resize(data.mz.size());
        for (size_t i = 0; i < data.mz.size(); ++i)
        {
          spectrum[i].setMZ(data.mz[i]);
          spectrum[i].setIntensity(static_cast<Peak1D::IntensityType>(data.intensities[i]));
        }

        exp.addSpectrum(std::move(spectrum));
      }

      endProgress();

      // Sort by RT
      exp.sortSpectra(true);
      exp.updateRanges();

      OPENMS_LOG_INFO << "ThermoRawFile: loaded " << exp.size() << " spectra\n";
    }
    catch (const openms::thermo_bridge::bridge_error& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        path, String("Thermo bridge error: ") + e.what());
    }
  }

} // namespace OpenMS

#endif // WITH_THERMO_RAW
