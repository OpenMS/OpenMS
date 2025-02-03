// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>

namespace OpenMS
{

  CachedmzML::CachedmzML() = default;

  CachedmzML::CachedmzML(const String& filename)
  {
    load_(filename);
  }

  CachedmzML::~CachedmzML()
  {
    ifs_.close();
  }

  CachedmzML::CachedmzML(const CachedmzML & rhs) :
    meta_ms_experiment_(rhs.meta_ms_experiment_),
    ifs_(rhs.filename_cached_.c_str(), std::ios::binary),
    filename_(rhs.filename_),
    spectra_index_(rhs.spectra_index_),
    chrom_index_(rhs.chrom_index_)
  {
  }

  void CachedmzML::load_(const String& filename)
  {
    filename_cached_ = filename + ".cached";
    filename_ = filename;

    // Create the index from the given file
    Internal::CachedMzMLHandler cache;
    cache.createMemdumpIndex(filename_cached_);
    spectra_index_ = cache.getSpectraIndex();
    chrom_index_ = cache.getChromatogramIndex();;

    // open the filestream
    ifs_.open(filename_cached_.c_str(), std::ios::binary);

    // load the meta data from disk
    FileHandler().loadExperiment(filename, meta_ms_experiment_, {OpenMS::FileTypes::MZML});
  }

  MSSpectrum CachedmzML::getSpectrum(Size id)
  {
    OPENMS_PRECONDITION(id < getNrSpectra(), "Id cannot be larger than number of spectra");

    if ( !ifs_.seekg(spectra_index_[id]) )
    {
      std::cerr << "Error while reading spectrum " << id << " - seekg created an error when trying to change position to " << spectra_index_[id] << "." << std::endl;
      std::cerr << "Maybe an invalid position was supplied to seekg, this can happen for example when reading large files (>2GB) on 32bit systems." << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Error while changing position of input stream pointer.", filename_cached_);
    }

    MSSpectrum s = meta_ms_experiment_.getSpectrum(id);
    Internal::CachedMzMLHandler::readSpectrum(s, ifs_);
    return s;
  }

  MSChromatogram CachedmzML::getChromatogram(Size id)
  {
    OPENMS_PRECONDITION(id < getNrChromatograms(), "Id cannot be larger than number of chromatograms");

    if ( !ifs_.seekg(chrom_index_[id]) )
    {
      std::cerr << "Error while reading chromatogram " << id << " - seekg created an error when trying to change position to " << chrom_index_[id] << "." << std::endl;
      std::cerr << "Maybe an invalid position was supplied to seekg, this can happen for example when reading large files (>2GB) on 32bit systems." << std::endl;
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Error while changing position of input stream pointer.", filename_cached_);
    }

    MSChromatogram c = meta_ms_experiment_.getChromatogram(id);
    Internal::CachedMzMLHandler::readChromatogram(c, ifs_);
    return c;
  }

  size_t CachedmzML::getNrSpectra() const
  {
    return meta_ms_experiment_.size();
  }

  size_t CachedmzML::getNrChromatograms() const
  {
    return meta_ms_experiment_.getChromatograms().size();
  }

  void CachedmzML::store(const String& filename, const PeakMap& map)
  {
    Internal::CachedMzMLHandler().writeMemdump(map, filename + ".cached");
    Internal::CachedMzMLHandler().writeMetadata_x(map, filename, true);
  }

  void CachedmzML::load(const String& filename, CachedmzML& map)
  {
    map.load_(filename);
  }

}

