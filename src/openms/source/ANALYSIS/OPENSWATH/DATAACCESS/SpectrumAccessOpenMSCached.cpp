// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>

#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>

namespace OpenMS
{

  SpectrumAccessOpenMSCached::SpectrumAccessOpenMSCached(const String& filename) :
    CachedmzML(filename)
  {
  }

  SpectrumAccessOpenMSCached::~SpectrumAccessOpenMSCached() = default;

  SpectrumAccessOpenMSCached::SpectrumAccessOpenMSCached(const SpectrumAccessOpenMSCached & rhs) :
    CachedmzML(rhs)
  {
    // this only copies the indices and meta-data
  }

  boost::shared_ptr<OpenSwath::ISpectrumAccess> SpectrumAccessOpenMSCached::lightClone() const
  {
    return boost::shared_ptr<SpectrumAccessOpenMSCached>(new SpectrumAccessOpenMSCached(*this));
  }

  OpenSwath::SpectrumPtr SpectrumAccessOpenMSCached::getSpectrumById(int id)
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");

    int ms_level = -1;
    double rt = -1.0;

    if ( !ifs_.seekg(spectra_index_[id]) )
    {
      std::cerr << "Error while reading spectrum " << id << " - seekg created an error when trying to change position to " << spectra_index_[id] << "." << std::endl;
      std::cerr << "Maybe an invalid position was supplied to seekg, this can happen for example when reading large files (>2GB) on 32bit systems." << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Error while changing position of input stream pointer.", filename_cached_);
    }

    OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
    sptr->getDataArrays() = Internal::CachedMzMLHandler::readSpectrumFast(ifs_, ms_level, rt);

    return sptr;
  }

  OpenSwath::SpectrumMeta SpectrumAccessOpenMSCached::getSpectrumMetaById(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrSpectra(), "Id cannot be larger than number of spectra");

    OpenSwath::SpectrumMeta meta;
    meta.RT = meta_ms_experiment_[id].getRT();
    meta.ms_level = meta_ms_experiment_[id].getMSLevel();
    return meta;
  }

  OpenSwath::ChromatogramPtr SpectrumAccessOpenMSCached::getChromatogramById(int id)
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of chromatograms");

    if ( !ifs_.seekg(chrom_index_[id]) )
    {
      std::cerr << "Error while reading chromatogram " << id << " - seekg created an error when trying to change position to " << chrom_index_[id] << "." << std::endl;
      std::cerr << "Maybe an invalid position was supplied to seekg, this can happen for example when reading large files (>2GB) on 32bit systems." << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Error while changing position of input stream pointer.", filename_cached_);
    }

    OpenSwath::ChromatogramPtr cptr(new OpenSwath::Chromatogram);
    cptr->getDataArrays() = Internal::CachedMzMLHandler::readChromatogramFast(ifs_);
    return cptr;
  }

  std::vector<std::size_t> SpectrumAccessOpenMSCached::getSpectraByRT(double RT, double deltaRT) const
  {
    OPENMS_PRECONDITION(deltaRT >= 0, "Delta RT needs to be a positive number");

    // we first perform a search for the spectrum that is past the
    // beginning of the RT domain. Then we add this spectrum and try to add
    // further spectra as long as they are below RT + deltaRT.
    std::vector<std::size_t> result;
    auto spectrum = meta_ms_experiment_.RTBegin(RT - deltaRT);
    if (spectrum == meta_ms_experiment_.end()) return result;

    result.push_back(std::distance(meta_ms_experiment_.begin(), spectrum));
    spectrum++;

    while (spectrum != meta_ms_experiment_.end() && spectrum->getRT() < RT + deltaRT)
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
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of spectra");
    return meta_ms_experiment_.getChromatograms()[id];
  }

  std::string SpectrumAccessOpenMSCached::getChromatogramNativeID(int id) const
  {
    OPENMS_PRECONDITION(id >= 0, "Id needs to be larger than zero");
    OPENMS_PRECONDITION(id < (int)getNrChromatograms(), "Id cannot be larger than number of spectra");
    return meta_ms_experiment_.getChromatograms()[id].getNativeID();
  }

} //end namespace OpenMS

