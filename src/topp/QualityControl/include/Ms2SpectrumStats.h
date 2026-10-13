// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>


namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  class TransformationDescription;

  /**
    @brief  QC metric to determine the number of MS2 scans per MS1 scan over RT

    Ms2SpectrumStats collects data from MS2 scans and stores the result at the identifications of the FeatureMap, or at
    new observations for the MS2 scans without identification.

    The following meta-values are computed:
    "ScanEventNumber": consecutive number of each MS2 scan after the preceding MS1 scan
    "identified": All identifications of the FeatureMap are marked with 1 and all unidentified MS2-Spectra with 0.
    "ion_injection_time": from MS2 spectrum
    "activation_method": from MS2 spectrum
    "total_ion_count": summed intensity from MS2 spectrum
    "base_peak_intensity": highest intensity from MS2 spectrum

    "FWHM": RT peak width for all assigned PIs (if provided)

    **/
  class Ms2SpectrumStats : public QCBase
  {
  public:
    struct ScanEvent {
      ScanEvent(UInt32 sem, bool ms2) : scan_event_number(sem), ms2_presence(ms2)
      {
      }
      UInt32 scan_event_number;
      bool ms2_presence;
    };

    /// Constructor
    Ms2SpectrumStats() = default;

    /// Destructor
    virtual ~Ms2SpectrumStats() = default;

    /**
      @brief Calculate the ScanEventNumber, find all unidentified MS2-Spectra and return them as observations,
             write meta values "ScanEventNumber" and "identified" in the identifications of the features and the unassigned ones
             (see QCBase::annotateIdentifications(); maps with peptide identifications are converted for this and back).
      @param[in] exp Imported calibrated MzML file as MSExperiment
      @param[in,out] features Imported featureXML file after FDR as FeatureMap
      @param[in] map_to_spectrum Map to find index of spectrum given by the spectrum reference of an identification
      @return observations of the unidentified MS2-Spectra (retention time, precursor m/z, spectrum reference and meta values),
              for identifications without matches
      @throws MissingInformation If exp is empty
      @throws InvalidParameter An identification has no spectrum reference
    **/
    std::vector<IdentificationData::Observation> compute(const MSExperiment& exp, FeatureMap& features, const QCBase::SpectraMap& map_to_spectrum);

    /// returns the name of the metric
    const std::string& getName() const override;
    /// define the required input file: featureXML after FDR (=POSTFDRFEAT), MzML-file (MSExperiment) with all MS2-Spectra (=RAWMZML)
    Status requirements() const override;

  private:
    /// name of the metric
    const std::string name_ = "Ms2SpectrumStats";

    /// ms2_included_ contains for every spectrum the information "ScanEventNumber" and presence MS2-scan in PeptideIDs
    std::vector<ScanEvent> ms2_included_ {};

    /// compute "ScanEventNumber" for every spectrum: MS1=0, MS2=1-n, write into ms2_included_
    void setScanEventNumber_(const MSExperiment& exp);

    /// set ms2_included_ bool to true, if PeptideID exist and set "ScanEventNumber" for every PeptideID
    void setPresenceAndScanEventNumber_(QCBase::AnnotatedIdentification& id, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum);

    /// return all unidentified MS2-Scans as observations, with their RT, precursor m/z, "ScanEventNumber" and spectrum meta values
    std::vector<IdentificationData::Observation> getUnidentifiedObservations_(const MSExperiment& exp);

    /// calculate highest intensity (base peak intensity)
    static MSSpectrum::PeakType::IntensityType getBPI_(const MSSpectrum& spec);
  };
} // namespace OpenMS
