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

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>

#include <fstream>

namespace OpenMS
{

  /**
    @brief An class that uses on-disk caching to read and write spectra and chromatograms

    This class provides functions to read and write spectra and chromatograms
    to disk using a time-efficient format. Reading the data items from disk can
    be very fast and done in random order (once the in-memory index is built
    for the file).

  */
  class OPENMS_DLLAPI CachedmzML
  {

public:

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    CachedmzML();

    CachedmzML(const String& filename);

    /// Copy constructor
    CachedmzML(const CachedmzML & rhs);

    /// Default destructor
    ~CachedmzML();
    //@}

    MSSpectrum getSpectrum(Size id);

    MSChromatogram getChromatogram(Size id);

    size_t getNrSpectra() const;

    size_t getNrChromatograms() const;

    const MSExperiment& getMetaData() const
    {
      return meta_ms_experiment_;
    }

    /**
      @brief Stores a map in a cached MzML file.

      @p filename The data location (ends in .mzML)
      @p map has to be an MSExperiment or have the same interface.

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    static void store(const String& filename, const PeakMap& map);

    /**
      @brief Loads a map from a cached MzML file

      @p filename The data location (ends in .mzML, expects an adjacent .mzML.cached file)
      @p map A CachedmzML result object

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    static void load(const String& filename, CachedmzML& map);

protected:

    void load_(const String& filename);

    /// Meta data
    MSExperiment meta_ms_experiment_;

    /// Internal filestream 
    std::ifstream ifs_;

    /// Name of the mzML file
    String filename_;

    /// Name of the cached mzML file
    String filename_cached_;

    /// Indices
    std::vector<std::streampos> spectra_index_;
    std::vector<std::streampos> chrom_index_;

  };
}

