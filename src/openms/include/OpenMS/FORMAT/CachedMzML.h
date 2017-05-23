// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_FORMAT_CACHEDMZML_H
#define OPENMS_FORMAT_CACHEDMZML_H

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <fstream>

#define CACHED_MZML_FILE_IDENTIFIER 8093

namespace OpenMS
{

  /**
    @brief An class that uses on-disk caching to read and write spectra and chromatograms

    This class provides functions to read and write spectra and chromatograms
    to disk using a time-efficient format. Reading the data items from disk can
    be very fast and done in random order (once the in-memory index is built
    for the file).

  */
  class OPENMS_DLLAPI CachedmzML :
    public ProgressLogger
  {
    int int_field_;
    double dbl_field_;

public:

    typedef PeakMap MapType;
    typedef MSSpectrum<Peak1D> SpectrumType;
    typedef MSChromatogram<ChromatogramPeak> ChromatogramType;

    // using double precision to store all data (has to agree with type of BinaryDataArrayPtr)
    typedef double DatumSingleton;

    typedef std::vector<DatumSingleton> Datavector;

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    CachedmzML();

    /// Default destructor
    ~CachedmzML();

    /// Assignment operator
    CachedmzML& operator=(const CachedmzML& rhs);
    //@}

    /** @name Read / Write a complete mass spectrometric experiment (or its meta data)
    */
    //@{

    /// Write complete spectra as a dump to the disk
    void writeMemdump(MapType& exp, String out);

    /// Write only the meta data of an MSExperiment
    void writeMetadata(MapType exp, String out_meta, bool addCacheMetaValue=false);

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

      @throws Exception::ParseError is thrown if the spectrum size cannot be read
    */
    static inline void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1,
                                        OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs, int& ms_level,
                                        double& rt)
    {
      Size spec_size = -1;
      ifs.read((char*) &spec_size, sizeof(spec_size));
      ifs.read((char*) &ms_level, sizeof(ms_level));
      ifs.read((char*) &rt, sizeof(rt));

      if ( static_cast<int>(spec_size) < 0)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Read an invalid spectrum length, something is wrong here. Aborting.", "filestream");
      }

      data1->data.resize(spec_size);
      data2->data.resize(spec_size);

      if (spec_size > 0)
      {
        ifs.read((char*) &(data1->data)[0], spec_size * sizeof(double));
        ifs.read((char*) &(data2->data)[0], spec_size * sizeof(double));
      }
    }

    /**
      @brief fast access to a chromatogram (a direct copy of the data into the provided arrays)

      @throws Exception::ParseError is thrown if the chromatogram size cannot be read
    */
    static inline void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1,
                                            OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs)
    {
      Size spec_size = -1;
      ifs.read((char*) &spec_size, sizeof(spec_size));

      if ( static_cast<int>(spec_size) < 0)
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Read an invalid chromatogram length, something is wrong here. Aborting.", "filestream");
      }

      data1->data.resize(spec_size);
      data2->data.resize(spec_size);
      ifs.read((char*) &(data1->data)[0], spec_size * sizeof(double));
      ifs.read((char*) &(data2->data)[0], spec_size * sizeof(double));
    }
    //@}

protected:

    /// read a single spectrum directly into a datavector (assuming file is already at the correct position)
    void readSpectrum_(Datavector& data1, Datavector& data2, std::ifstream& ifs, int& ms_level, double& rt) const;

    /// read a single chromatogram directly into a datavector (assuming file is already at the correct position)
    void readChromatogram_(Datavector& data1, Datavector& data2, std::ifstream& ifs) const;

    /// read a single spectrum directly into an OpenMS MSSpectrum (assuming file is already at the correct position)
    void readSpectrum_(SpectrumType& spectrum, std::ifstream& ifs) const;

    /// read a single chromatogram directly into an OpenMS MSChromatograms (assuming file is already at the correct position)
    void readChromatogram_(ChromatogramType& chromatogram, std::ifstream& ifs) const;

    /// write a single spectrum to filestream
    void writeSpectrum_(const SpectrumType& spectrum, std::ofstream& ofs);

    /// write a single chromatogram to filestream
    void writeChromatogram_(const ChromatogramType& chromatogram, std::ofstream& ofs);

    /// Members
    std::vector<std::streampos> spectra_index_;
    std::vector<std::streampos> chrom_index_;

  };
}
#endif // OPENMS_FORMAT_CACHEDMZML_H

