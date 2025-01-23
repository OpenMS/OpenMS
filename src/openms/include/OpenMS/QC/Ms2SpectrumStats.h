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


namespace OpenMS
{
  class FeatureMap;
  class MSExperiment;
  class PeptideIdentification;
  class TransformationDescription;

  /**
    @brief  QC metric to determine the number of MS2 scans per MS1 scan over RT

    Ms2SpectrumStats collects data from MS2 scans and stores the result into PeptideIdentifications,
    which already exist in the FeatureMap, or are newly created as empty PeptideIdentifications (with no sequence).

    The following meta-values are computed:
    "ScanEventNumber": consecutive number of each MS2 scan after the preceding MS1 scan
    "identified": All PeptideIdentifications of the FeatureMap are marked with '+' and all unidentified MS2-Spectra with '-'.
    "ion_injection_time": from MS2 spectrum
    "activation_method": from MS2 spectrum
    "total_ion_count": summed intensity from MS2 spectrum
    "base_peak_intensity": highest intensity from MS2 spectrum

    "FWHM": RT peak width for all assigned PIs (if provided)

    **/
  class OPENMS_DLLAPI Ms2SpectrumStats : public QCBase
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
      @brief Calculate the ScanEventNumber, find all unidentified MS2-Spectra and add them to unassigned PeptideIdentifications,
             write meta values "ScanEventNumber" and "identified" in PeptideIdentification.
      @param exp Imported calibrated MzML file as MSExperiment
      @param features Imported featureXML file after FDR as FeatureMap
      @param map_to_spectrum Map to find index of spectrum given by meta value at PepID
      @return unassigned peptide identifications newly generated from unidentified MS2-Spectra
      @throws MissingInformation If exp is empty
      @throws InvalidParameter PeptideID is missing meta value 'spectrum_reference'
    **/
    std::vector<PeptideIdentification> compute(const MSExperiment& exp, FeatureMap& features, const QCBase::SpectraMap& map_to_spectrum);

    /// returns the name of the metric
    const String& getName() const override;
    /// define the required input file: featureXML after FDR (=POSTFDRFEAT), MzML-file (MSExperiment) with all MS2-Spectra (=RAWMZML)
    Status requirements() const override;

  private:
    /// name of the metric
    const String name_ = "Ms2SpectrumStats";

    /// ms2_included_ contains for every spectrum the information "ScanEventNumber" and presence MS2-scan in PeptideIDs
    std::vector<ScanEvent> ms2_included_ {};

    /// compute "ScanEventNumber" for every spectrum: MS1=0, MS2=1-n, write into ms2_included_
    void setScanEventNumber_(const MSExperiment& exp);

    /// set ms2_included_ bool to true, if PeptideID exist and set "ScanEventNumber" for every PeptideID
    void setPresenceAndScanEventNumber_(PeptideIdentification& peptide_ID, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum);

    /// return all unidentified MS2-Scans as unassignedPeptideIDs, these contain only Information about RT and "ScanEventNumber"
    std::vector<PeptideIdentification> getUnassignedPeptideIdentifications_(const MSExperiment& exp);

    /// calculate highest intensity (base peak intensity)
    static MSSpectrum::PeakType::IntensityType getBPI_(const MSSpectrum& spec);
  };
} // namespace OpenMS
