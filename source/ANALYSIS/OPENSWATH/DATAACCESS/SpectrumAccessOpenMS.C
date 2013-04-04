// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>

namespace OpenMS
{
  SpectrumAccessOpenMS::SpectrumAccessOpenMS(MSExperimentType& ms_experiment)
  {
    // store raw ptr to the experiment
    ms_experiment_ = &ms_experiment;
  }

  SpectrumAccessOpenMS::~SpectrumAccessOpenMS()
  {
  }

  OpenSwath::SpectrumPtr SpectrumAccessOpenMS::getSpectrumById(int id) const
  {
    const MSSpectrumType& spectrum = (*ms_experiment_)[id];
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    for (MSSpectrumType::const_iterator it = spectrum.begin(); it != spectrum.end(); it++)
    {
      mz_array->data.push_back(it->getMZ());
      intensity_array->data.push_back(it->getIntensity());
    }

    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    sptr->setMZArray(mz_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  OpenSwath::SpectrumMeta SpectrumAccessOpenMS::getSpectrumMetaById(int id) const
  {
    OpenSwath::SpectrumMeta meta;
    meta.RT = (*ms_experiment_)[id].getRT();
    meta.ms_level = (*ms_experiment_)[id].getMSLevel();
    return meta;
  }

  OpenSwath::ChromatogramPtr SpectrumAccessOpenMS::getChromatogramById(int id) const
  {
    const MSChromatogramType& chromatogram = ms_experiment_->getChromatograms()[id];
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr rt_array(new OpenSwath::BinaryDataArray);
    for (MSChromatogramType::const_iterator it = chromatogram.begin(); it != chromatogram.end(); it++)
    {
      rt_array->data.push_back(it->getRT());
      intensity_array->data.push_back(it->getIntensity());
    }

    // push back rt first, then intensity.
    // FEATURE (hroest) annotate which is which
    std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
    binaryDataArrayPtrs.push_back(rt_array);
    binaryDataArrayPtrs.push_back(intensity_array);

    OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram);
    cptr->binaryDataArrayPtrs = binaryDataArrayPtrs;
    return cptr;
  }

  std::vector<std::size_t> SpectrumAccessOpenMS::getSpectraByRT(double RT, double deltaRT) const
  {
    OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");

    // we first perform a search for the spectrum that is past the
    // beginning of the RT domain. Then we add this spectrum and try to add
    // further spectra as long as they are below RT + deltaRT.
    MSExperimentType::Iterator spectrum = ms_experiment_->RTBegin(RT - deltaRT);
    std::vector<std::size_t> result;
    result.push_back(std::distance(ms_experiment_->begin(), spectrum));
    spectrum++;
    while (spectrum->getRT() <= RT + deltaRT && spectrum != ms_experiment_->end())
    {
      result.push_back(spectrum - ms_experiment_->begin());
      spectrum++;
    }
    return result;
  }

  size_t SpectrumAccessOpenMS::getNrChromatograms() const
  {
    return ms_experiment_->getChromatograms().size();
  }

  ChromatogramSettings SpectrumAccessOpenMS::getChromatogramMetaInfo(int id) const
  {
    return ms_experiment_->getChromatograms()[id];
  }

  std::string SpectrumAccessOpenMS::getChromatogramNativeID(int id) const
  {
    return ms_experiment_->getChromatograms()[id].getNativeID();
  }

  size_t SpectrumAccessOpenMS::getNrSpectra() const
  {
    return ms_experiment_->size();
  }

  SpectrumSettings SpectrumAccessOpenMS::getSpectraMetaInfo(int id) const
  {
    return (*ms_experiment_)[id];
  }

} //end namespace OpenMS
