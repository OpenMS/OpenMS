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

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <fstream>

#define CACHED_MZML_FILE_IDENTIFIER 8094

namespace OpenMS
{

namespace Internal
{

  /**
    @brief An class that uses on-disk caching to read and write spectra and chromatograms

    This class provides functions to read and write spectra and chromatograms
    to disk using a time-efficient format. Reading the data items from disk can
    be very fast and done in random order (once the in-memory index is built
    for the file).

  */
  class OPENMS_DLLAPI CachedMzMLHandler :
    public ProgressLogger
  {
    typedef int IntType;
    typedef double DoubleType;

public:

    typedef PeakMap MapType;
    typedef MSSpectrum SpectrumType;
    typedef MSChromatogram ChromatogramType;

    // using double precision to store all data (has to agree with type of BinaryDataArrayPtr)
    typedef double DatumSingleton;

    typedef std::vector<DatumSingleton> Datavector;

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    CachedMzMLHandler();

    /// Default destructor
    ~CachedMzMLHandler();

    /// Assignment operator
    CachedMzMLHandler& operator=(const CachedMzMLHandler& rhs);
    //@}

    /** @name Read / Write a complete mass spectrometric experiment (or its meta data)
    */
    //@{

    /// Write complete spectra as a dump to the disk
    void writeMemdump(const MapType& exp, const String& out) const;

    /// Write only the meta data of an MSExperiment
    void writeMetadata(MapType exp, String out_meta, bool addCacheMetaValue=false);

    /// Write only the meta data of an MSExperiment
    void writeMetadata_x(const MapType& exp, const String& out_meta, bool addCacheMetaValue=false);

    /// Read all spectra from a dump from the disk
    void readMemdump(MapType& exp_reading, String filename) const;
    //@}

    /** @name Access and creation of the binary indices
    */
    //@{
    /// Create an index on the location of all the spectra and chromatograms
    void createMemdumpIndex(String filename);

    /// Access to a constant copy of the binary spectra index
    const std::vector<std::streampos>& getSpectraIndex() const;

    /// Access to a constant copy of the binary chromatogram index
    const std::vector<std::streampos>& getChromatogramIndex() const;
    //@}

    /** @name Direct access to a single Spectrum or Chromatogram
    */
    //@{

    /**
      @brief fast access to a spectrum (a direct copy of the data into the provided arrays)

      @param data1 First data array (m/z)
      @param data2 Second data array (Intensity)
      @param ms_level Output parameter to store the MS level of the spectrum (1, 2, 3 ...)
      @param rt Output parameter to store the retention time of the spectrum

      @throws Exception::ParseError is thrown if the spectrum cannot be read
    */
    static inline void readSpectrumFast(OpenSwath::BinaryDataArrayPtr& data1,
                                        OpenSwath::BinaryDataArrayPtr& data2,
                                        std::ifstream& ifs, 
                                        int& ms_level,
                                        double& rt)
    {
      std::vector<OpenSwath::BinaryDataArrayPtr> data = readSpectrumFast(ifs, ms_level, rt);
      data1 = data[0];
      data2 = data[1];
    }

    /**
      @brief Fast access to a spectrum

      @param ifs Input file stream (moved to the correct position)
      @param ms_level Output parameter to store the MS level of the spectrum (1, 2, 3 ...)
      @param rt Output parameter to store the retention time of the spectrum

      @throws Exception::ParseError is thrown if the spectrum cannot be read
    */
    static std::vector<OpenSwath::BinaryDataArrayPtr> readSpectrumFast(std::ifstream& ifs, int& ms_level, double& rt);

    /**
      @brief Fast access to a chromatogram

      @param data1 First data array (RT)
      @param data2 Second data array (Intensity)

      @throws Exception::ParseError is thrown if the chromatogram size cannot be read
    */
    static inline void readChromatogramFast(OpenSwath::BinaryDataArrayPtr& data1,
                                            OpenSwath::BinaryDataArrayPtr& data2, std::ifstream& ifs)
    {
      std::vector<OpenSwath::BinaryDataArrayPtr> data = readChromatogramFast(ifs);
      data1 = data[0];
      data2 = data[1];
    }

    /**
      @brief Fast access to a chromatogram

      @param ifs Input file stream (moved to the correct position)

      @throws Exception::ParseError is thrown if the chromatogram size cannot be read
    */
    static std::vector<OpenSwath::BinaryDataArrayPtr> readChromatogramFast(std::ifstream& ifs);
    //@}

    /**
      @brief Read a single spectrum directly into an OpenMS MSSpectrum (assuming file is already at the correct position)

      @param spectrum Output spectrum
      @param ifs Input file stream (moved to the correct position)

      @throws Exception::ParseError is thrown if the chromatogram size cannot be read
    */
    static void readSpectrum(SpectrumType& spectrum, std::ifstream& ifs);

    /**
      @brief Read a single chromatogram directly into an OpenMS MSChromatogram (assuming file is already at the correct position)

      @param chromatogram Output chromatogram
      @param ifs Input file stream (moved to the correct position)

      @throws Exception::ParseError is thrown if the chromatogram size cannot be read
    */
    static void readChromatogram(ChromatogramType& chromatogram, std::ifstream& ifs);

protected:

    /// write a single spectrum to filestream
    void writeSpectrum_(const SpectrumType& spectrum, std::ofstream& ofs) const;

    /// write a single chromatogram to filestream
    void writeChromatogram_(const ChromatogramType& chromatogram, std::ofstream& ofs) const;

    /// helper method for fast reading of spectra and chromatograms
    static inline void readDataFast_(std::ifstream& ifs, std::vector<OpenSwath::BinaryDataArrayPtr>& data, const Size& data_size, 
      const Size& nr_float_arrays);

    /// Members
    std::vector<std::streampos> spectra_index_;
    std::vector<std::streampos> chrom_index_;

  };
}
}

