// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/Ms2SpectrumStats.h>
#include <OpenMS/QC/QCBase.h>

using namespace std;

namespace OpenMS
{
  // check which MS2-Spectra of a mzml-file (MSExperiment) are identified (and therefor have a entry in the featureMap)
  // MS2 spectra without mate are returned as vector of unassigned PeptideIdentifications (with empty sequence but some metavalue)
  std::vector<PeptideIdentification> Ms2SpectrumStats::compute(const MSExperiment& exp, FeatureMap& features, const QCBase::SpectraMap& map_to_spectrum)
  {
    if (exp.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The mzml file / MSExperiment must not be empty.\n");
    }

    setScanEventNumber_(exp);
    // if MS2-spectrum PeptideIdentifications found ->  ms2_included_ nullptr to PepID pointer
    std::function<void(PeptideIdentification&)> l_f = [&exp, this, &map_to_spectrum](PeptideIdentification& pep_id) { setPresenceAndScanEventNumber_(pep_id, exp, map_to_spectrum); };
    features.applyFunctionOnPeptideIDs(l_f);

    // if Ms2-spectrum not identified, add to unassigned PeptideIdentification without ID, contains only RT, mz and some meta values
    return getUnassignedPeptideIdentifications_(exp);
  }


  void Ms2SpectrumStats::setScanEventNumber_(const MSExperiment& exp)
  {
    ms2_included_.clear();
    ms2_included_.reserve(exp.size());
    UInt32 scan_event_number {0};
    for (const MSSpectrum& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() == 1)
      { // reset
        scan_event_number = 0;
        ms2_included_.emplace_back(scan_event_number, false);
      }
      else if (spec.getMSLevel() == 2)
      {
        ++scan_event_number;
        ms2_included_.emplace_back(scan_event_number, false);
      }
    }
  }

  void annotatePepIDfromSpectrum_(const MSSpectrum& spectrum, PeptideIdentification& peptide_ID)
  {
    if (!spectrum.getAcquisitionInfo().empty() && spectrum.getAcquisitionInfo()[0].metaValueExists("MS:1000927"))
    {
      peptide_ID.setMetaValue("ion_injection_time", spectrum.getAcquisitionInfo()[0].getMetaValue("MS:1000927"));
    }
    if (!spectrum.getPrecursors().empty() && !spectrum.getPrecursors()[0].getActivationMethods().empty())
    {
      peptide_ID.setMetaValue("activation_method", Precursor::NamesOfActivationMethodShort[*spectrum.getPrecursors()[0].getActivationMethods().begin()]);
    }
  }

  // marks all seen (unassigned-)PeptideIdentifications in vector ms2_included
  void Ms2SpectrumStats::setPresenceAndScanEventNumber_(PeptideIdentification& peptide_ID, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum)
  {
    if (!peptide_ID.metaValueExists("spectrum_reference"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum reference annotated at peptide identification!");
    }

    UInt64 index = map_to_spectrum.at(peptide_ID.getSpectrumReference());
    const MSSpectrum& spectrum = exp[index];

    if (spectrum.getMSLevel() == 2)
    {
      ms2_included_[index].ms2_presence = true;
      peptide_ID.setMetaValue("ScanEventNumber", ms2_included_[index].scan_event_number);
      peptide_ID.setMetaValue("identified", 1);
      peptide_ID.setMetaValue("total_ion_count", spectrum.calculateTIC());
      peptide_ID.setMetaValue("base_peak_intensity", getBPI_(spectrum));
      annotatePepIDfromSpectrum_(spectrum, peptide_ID); // ion_injection_time and activation_method
    }
  }

  std::vector<PeptideIdentification> Ms2SpectrumStats::getUnassignedPeptideIdentifications_(const MSExperiment& exp)
  {
    std::vector<PeptideIdentification> result;
    for (auto it = ms2_included_.begin(); it != ms2_included_.end(); ++it)
    {
      if (it->ms2_presence)
      {
        continue;
      }
      const MSSpectrum& spec = exp.getSpectra()[distance(ms2_included_.begin(), it)];
      if (spec.getMSLevel() != 2)
      {
        continue;
      }
      PeptideIdentification unidentified_MS2;
      unidentified_MS2.setRT(spec.getRT());
      unidentified_MS2.setMetaValue("ScanEventNumber", (*it).scan_event_number);
      unidentified_MS2.setMetaValue("identified", 0);
      unidentified_MS2.setMZ(spec.getPrecursors()[0].getMZ());
      unidentified_MS2.setMetaValue("total_ion_count", spec.calculateTIC());
      unidentified_MS2.setMetaValue("base_peak_intensity", getBPI_(spec));
      unidentified_MS2.setSpectrumReference( spec.getNativeID());
      annotatePepIDfromSpectrum_(spec, unidentified_MS2); // ion_injection_time and activation_method
      result.push_back(unidentified_MS2);
    }
    return result;
  }

  // calculate maximal and summed intensity
  MSSpectrum::PeakType::IntensityType Ms2SpectrumStats::getBPI_(const MSSpectrum& spec)
  {
    PeakSpectrum::PeakType::IntensityType bpi {0};
    auto it = spec.getBasePeak();
    if (it != spec.end())
    {
      bpi = it->getIntensity();
    }
    return bpi;
  }

  // returns the name of the metric
  const String& Ms2SpectrumStats::getName() const
  {
    return name_;
  }

  // required input files
  QCBase::Status Ms2SpectrumStats::requirements() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }

} // namespace OpenMS
