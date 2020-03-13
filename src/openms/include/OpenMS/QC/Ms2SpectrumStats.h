// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/QC/QCBase.h>

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>


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
    struct ScanEvent
    {
      ScanEvent(UInt32 sem, bool ms2)
        : scan_event_number(sem), ms2_presence(ms2) 
      {}
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
    Status requires() const override;

  private:
    /// name of the metric
    const String name_ = "Ms2SpectrumStats";
    
    /// ms2_included_ contains for every spectrum the information "ScanEventNumber" and presence MS2-scan in PeptideIDs
    std::vector<ScanEvent> ms2_included_{};

    /// compute "ScanEventNumber" for every spectrum: MS1=0, MS2=1-n, write into ms2_included_
    void setScanEventNumber_(const MSExperiment& exp);

    /// set ms2_included_ bool to true, if PeptideID exist and set "ScanEventNumber" for every PeptideID
    void setPresenceAndScanEventNumber_(PeptideIdentification& peptide_ID, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum);

    /// return all unidentified MS2-Scans as unassignedPeptideIDs, these contain only Information about RT and "ScanEventNumber"
    std::vector<PeptideIdentification> getUnassignedPeptideIdentifications_(const MSExperiment& exp);

    /// calculate highest intensity (base peak intensity)
    static MSSpectrum::PeakType::IntensityType getBPI_(const MSSpectrum& spec);
  };
}
