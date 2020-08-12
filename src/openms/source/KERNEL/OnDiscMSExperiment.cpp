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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <OpenMS/FORMAT/MzMLFile.h>

namespace OpenMS
{

  void OnDiscMSExperiment::loadMetaData_(const String& filename)
  {
    meta_ms_experiment_ = boost::shared_ptr< PeakMap >(new PeakMap);

    MzMLFile f;
    PeakFileOptions options = f.getOptions();
    options.setFillData(false);
    f.setOptions(options);
    f.load(filename, *meta_ms_experiment_.get());
  }

  MSChromatogram OnDiscMSExperiment::getMetaChromatogramById_(const std::string& id)
  {
    if (chromatograms_native_ids_.empty())
    {
      for (Size k = 0; k < meta_ms_experiment_->getChromatograms().size(); k++)
      {
        chromatograms_native_ids_.emplace(meta_ms_experiment_->getChromatograms()[k].getNativeID(), k);
      }
    }

    if (chromatograms_native_ids_.find(id) == chromatograms_native_ids_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("Could not find chromatogram with id '") + id + "'.");
    }
    return meta_ms_experiment_->getChromatogram(chromatograms_native_ids_[id]);
  }

  MSChromatogram OnDiscMSExperiment::getChromatogramByNativeId(const std::string& id)
  {
    if (!meta_ms_experiment_)
    {
      MSChromatogram chromatogram;
      indexed_mzml_file_.getMSChromatogramByNativeId(id, chromatogram);
      return chromatogram;
    }

    MSChromatogram chromatogram = getMetaChromatogramById_(id);
    indexed_mzml_file_.getMSChromatogramByNativeId(id, chromatogram);
    return chromatogram;
  }

  MSSpectrum OnDiscMSExperiment::getMetaSpectrumById_(const std::string& id)
  {
    if (spectra_native_ids_.empty())
    {
      for (Size k = 0; k < meta_ms_experiment_->getSpectra().size(); k++)
      {
        spectra_native_ids_.emplace(meta_ms_experiment_->getSpectra()[k].getNativeID(), k);
      }
    }

    if (spectra_native_ids_.find(id) == spectra_native_ids_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("Could not find spectrum with id '") + id + "'.");
    }
    return meta_ms_experiment_->getSpectrum(spectra_native_ids_[id]);
  }

  MSSpectrum OnDiscMSExperiment::getSpectrumByNativeId(const std::string& id)
  {
    if (!meta_ms_experiment_)
    {
      MSSpectrum spec;
      indexed_mzml_file_.getMSSpectrumByNativeId(id, spec);
      return spec;
    }

    MSSpectrum spec = getMetaSpectrumById_(id);
    indexed_mzml_file_.getMSSpectrumByNativeId(id, spec);
    return spec;
  }

} //namespace OpenMS

