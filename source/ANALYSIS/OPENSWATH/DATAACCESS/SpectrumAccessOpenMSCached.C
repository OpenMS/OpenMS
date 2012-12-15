// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>

namespace OpenMS
{
  OpenSwath::SpectrumPtr SpectrumAccessOpenMSCached::getSpectrumById(int id) const
  {
    if (cache_.getSpectraIndex().empty())
    {
      // remove const from the cache since we need to recalculate the index
      // and re-read the data.
      (const_cast<CachedmzML*>(&cache_))->createMemdumpIndex(filename_cached_);
    }
    OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    int ms_level = -1;
    double rt = -1.0;
    // FEATURE check if we can keep the filestream open -> risky if someone else
    // accesses the file in the meantime
    std::ifstream ifs_((filename_cached_).c_str(), std::ios::binary);
    ifs_.seekg(cache_.getSpectraIndex()[id]);
    cache_.readSpectrumFast(mz_array, intensity_array, ifs_, ms_level, rt);

    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    sptr->setMZArray(mz_array);
    sptr->setIntensityArray(intensity_array);
    return sptr;
  }

  OpenSwath::SpectrumMeta SpectrumAccessOpenMSCached::getSpectrumMetaById(int id) const
  {
    OpenSwath::SpectrumMeta meta;
    meta.RT = meta_ms_experiment_[id].getRT();
    meta.ms_level = meta_ms_experiment_[id].getMSLevel();
    return meta;
  }

  OpenSwath::ChromatogramPtr SpectrumAccessOpenMSCached::getChromatogramById(int id) const
  {
    if (cache_.getChromatogramIndex().empty())
    {
      // remove const from the cache since we need to recalculate the index
      // and re-read the data.
      (const_cast<CachedmzML*>(&cache_))->createMemdumpIndex(filename_cached_);
    }
    OpenSwath::BinaryDataArrayPtr rt_array(new OpenSwath::BinaryDataArray);
    OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
    std::ifstream ifs_((filename_cached_).c_str(), std::ios::binary);
    ifs_.seekg(cache_.getChromatogramIndex()[id]);
    cache_.readChromatogramFast(rt_array, intensity_array, ifs_);

    // push back rt first, then intensity.
    // FEATURE (hroest) annotate which is which
    std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
    binaryDataArrayPtrs.push_back(rt_array);
    binaryDataArrayPtrs.push_back(intensity_array);

    OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram);
    cptr->binaryDataArrayPtrs = binaryDataArrayPtrs;
    return cptr;
  }

  SpectrumAccessOpenMSCached::SpectrumAccessOpenMSCached(String filename)
  {
    filename_cached_ = filename + ".cached";
    // currently we re-open the filestream with each read access
    // std::ifstream ifs_((filename_cached_).c_str(), std::ios::binary);
    MzMLFile f;
    f.load(filename, meta_ms_experiment_);
    filename_ = filename;
  }

  SpectrumAccessOpenMSCached::~SpectrumAccessOpenMSCached()
  {
  }

  std::vector<std::size_t> SpectrumAccessOpenMSCached::getSpectraByRT(double RT, double deltaRT) const
  {
    OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");

    // we first perform a search for the spectrum that is past the
    // beginning of the RT domain. Then we add this spectrum and try to add
    // further spectra as long as they are below RT + deltaRT.
    std::vector<std::size_t> result;
    MSExperimentType::ConstIterator spectrum = meta_ms_experiment_.RTBegin(RT - deltaRT);
    result.push_back(std::distance(meta_ms_experiment_.begin(), spectrum));
    spectrum++;
    while (spectrum->getRT() < RT + deltaRT && spectrum != meta_ms_experiment_.end())
    {
      result.push_back(spectrum - meta_ms_experiment_.begin());
      spectrum++;
    }
    return result;
  }

  size_t SpectrumAccessOpenMSCached::getNrSpectra() const
  {
    return meta_ms_experiment_.size();
  }

  SpectrumSettings SpectrumAccessOpenMSCached::getSpectraMetaInfo(int id) const
  {
    return meta_ms_experiment_[id];
  }

  size_t SpectrumAccessOpenMSCached::getNrChromatograms() const
  {
    return meta_ms_experiment_.getChromatograms().size();
  }

  ChromatogramSettings SpectrumAccessOpenMSCached::getChromatogramMetaInfo(int id) const
  {
    return meta_ms_experiment_.getChromatograms()[id];
  }

  std::string SpectrumAccessOpenMSCached::getChromatogramNativeID(int id) const
  {
    return meta_ms_experiment_.getChromatograms()[id].getNativeID();
  }

} //end namespace OpenMS
