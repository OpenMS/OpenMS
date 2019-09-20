//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#include <OpenMS/QC/TopNoverRT.h>

#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{
  // check which MS2-Spectra of a mzml-file (MSExperiment) are identified (and therefor have a entry in the featureMap)
  // MS2 spectra without mate are returned as vector of unassignedPeptideIdentifications (with empty sequence but some metavalue)
  std::vector<PeptideIdentification> TopNoverRT::compute(const MSExperiment& exp, FeatureMap& features, const QCBase::SpectraMap& map_to_spectrum)
  {
    if (exp.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The mzml file / MSExperiment is empty.\n");
    }

    setScanEventNumber_(exp);
    // if MS2-spectrum PeptideIdentifications found ->  ms2_included_ nullptr to PepID pointer
    auto l_f = [&exp,this,&map_to_spectrum] (PeptideIdentification& pep_id)
    {
      setPresenceAndScanEventNumber_(pep_id, exp, map_to_spectrum);
    };
    iterateFeatureMap(features, l_f);
    for (auto& f : features)
    {
      if (f.metaValueExists("FWHM")) // from FF-Centroided
      {
        for (auto& pi : f.getPeptideIdentifications())
        {
          pi.setMetaValue("FWHM", f.getMetaValue("FWHM"));
        }
      }
      else if (f.metaValueExists("model_FWHM")) // from FF-Identification
      {
        for (auto& pi : f.getPeptideIdentifications())
        {
          pi.setMetaValue("FWHM", f.getMetaValue("model_FWHM")); // use 'FWHM' as target to make meta value unique for downstream processing
        }
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Metavalue 'FWHM' or 'model_FWHM' is missing for a feature in a FeatureMap. Please check your FeatureFinder reports FWHM using these metavalues or add a new mapping here.");
      }
    }

    // if Ms2-spectrum not identified, add to unassigned PeptideIdentification without ID, contains only RT and ScanEventNumber
    return getUnassignedPeptideIdentifications_(exp);
  }
  // if ms2 spetrum not included, add to unassignedPeptideIdentification, set m/z and RT values
  void TopNoverRT::setScanEventNumber_(const MSExperiment& exp)
  {
    ms2_included_.clear();
    ms2_included_.reserve(exp.size());
    UInt32 scan_event_number{ 0 };
    for (const MSSpectrum& spec : exp.getSpectra())
    {
      if (spec.getMSLevel() == 1)
      {
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
  void TopNoverRT::setPresenceAndScanEventNumber_(PeptideIdentification& peptide_ID, const MSExperiment& exp, const QCBase::SpectraMap& map_to_spectrum)
  {
    if (!peptide_ID.metaValueExists("spectrum_reference"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No spectrum reference annotated at peptide identification!");
    }
    
    UInt64 index = map_to_spectrum.at(peptide_ID.getMetaValue("spectrum_reference").toString());
    const MSSpectrum& spectrum = exp[index];

    if (spectrum.getMSLevel() == 2)
    {
      Peak1D::IntensityType bpi;
      Peak1D::IntensityType tic;
      getBPIandCIC_(spectrum, bpi, tic);
      ms2_included_[index].ms2_presence = true;
      peptide_ID.setMetaValue("ScanEventNumber", ms2_included_[index].scan_event_number);
      peptide_ID.setMetaValue("identified", 1);
      peptide_ID.setMetaValue("total_ion_count", tic);
      peptide_ID.setMetaValue("base_peak_intensity", bpi);
      annotatePepIDfromSpectrum_(spectrum, peptide_ID);
    }
  }

  std::vector<PeptideIdentification> TopNoverRT::getUnassignedPeptideIdentifications_(const MSExperiment& exp)
  {
    std::vector<PeptideIdentification> result;
    for (auto it = ms2_included_.begin(); it != ms2_included_.end(); ++it)
    {
      if (it->ms2_presence) continue;
      
      Size pos = distance(ms2_included_.begin(), it);
      if (exp[pos].getMSLevel() != 2) continue;

      PeptideIdentification unidentified_MS2;
      Peak1D::IntensityType bpi;
      Peak1D::IntensityType tic;
      getBPIandCIC_(exp.getSpectra()[pos], bpi, tic);
      unidentified_MS2.setRT(exp.getSpectra()[pos].getRT());
      unidentified_MS2.setMetaValue("ScanEventNumber", (*it).scan_event_number);
      unidentified_MS2.setMetaValue("identified", 0);
      unidentified_MS2.setMZ(exp.getSpectra()[pos].getPrecursors()[0].getMZ());
      unidentified_MS2.setMetaValue("total_ion_count", tic);
      unidentified_MS2.setMetaValue("base_peak_intensity", bpi);
      unidentified_MS2.setMetaValue("spectrum_reference", exp.getSpectra()[pos].getNativeID());
      annotatePepIDfromSpectrum_(exp.getSpectra()[pos], unidentified_MS2);
      result.push_back(unidentified_MS2);
    }
    return result;
  }

  // returns the name of the metric
  const String& TopNoverRT::getName() const
  {
    return name_;
  }

  // required input files
  QCBase::Status TopNoverRT::requires() const
  {
    return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
  }

  // calculate maximal and summed intensity
  void TopNoverRT::getBPIandCIC_(MSSpectrum const &spec,
                                  Peak1D::IntensityType& bpi,
                                  Peak1D::IntensityType& tic)
  {
    Peak1D::IntensityType peak_max{0};
    Peak1D::IntensityType sum{0};
    for (const Peak1D& peak : spec)
    {
      sum += peak.getIntensity();
      if (peak.getIntensity() > peak_max)
      {
        peak_max = peak.getIntensity();
      }
    }
    bpi = peak_max;
    tic = sum;
  }
}
