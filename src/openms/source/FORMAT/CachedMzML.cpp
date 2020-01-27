// --------------------------------------------------------------------------
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CachedMzML.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/FORMAT/HANDLERS/CachedMzMLHandler.h>

namespace OpenMS
{

  CachedmzML::CachedmzML()
  {
  }

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
    MzMLFile().load(filename, meta_ms_experiment_);
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

