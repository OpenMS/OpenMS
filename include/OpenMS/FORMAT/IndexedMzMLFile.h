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

#ifndef OPENMS_FORMAT_INDEXEDMZMLFILE_H
#define OPENMS_FORMAT_INDEXEDMZMLFILE_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/INTERFACES/ISpectrumAccess.h>

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>

#include <string>
#include <fstream>

//#define DEBUG_READER

namespace OpenMS
{

  /**
    @brief A class to read an indexedmzML file.

    This class implements access to an indexedmzML file and the contained spectra
    and chromatogram data through the getSpectrumById and getChromatogramById
    functions. It thus allows random access to spectra and chromatograms data
    without having to read the whole file into memory.

    Internally it uses the IndexedMzMLDecoder for initial parsing and
    extracting all the offsets of the <chromatogram> and <spectrum> tags. These
    offsets are stored as members of this class as well as the offset to the <indexList> element
  */
  class OPENMS_DLLAPI IndexedMzMLFile
  {
      /// Name of the file
      String filename_;
      /// Binary offsets to all spectra
      std::vector< std::pair<std::string, long> > spectra_offsets;
      /// Binary offsets to all chromatograms
      std::vector< std::pair<std::string, long> > chromatograms_offsets;
      /// offset to the <indexList> element
      long index_offset_;
      /// Whether spectra are written before chromatograms in this file
      bool spectra_before_chroms_;
      /// The current filestream (is opened upon construction)
      std::ifstream filestream; 
      /// Whether parsing the indexedmzML file was successful
      bool parsing_success_;

    /**
      @brief Try to parse the footer of the indexedmzML

      Upon success, the chromatogram and spectra offsets will be populated and
      parsing_success_ will be set to true.
    */
    void parseFooter(String filename);

    public:

    /**
      @brief Constructor

      Tries to parse the file, success can be checked with getParsingSuccess()
    */
    IndexedMzMLFile(String filename);

    /// Copy constructor
    IndexedMzMLFile(const IndexedMzMLFile & source);

    /// Destructor
    ~IndexedMzMLFile();

    /// Returns whether parsing was successful
    bool getParsingSuccess() const;

    /// Returns the number of spectra available
    size_t getNrSpectra() const;

    /// Returns the number of chromatograms available
    size_t getNrChromatograms() const;

    /// Returns the raw data for the spectrum at position "id"
    OpenMS::Interfaces::SpectrumPtr getSpectrumById(int id);

    /// Returns the raw data for the chromatogram at position "id"
    OpenMS::Interfaces::ChromatogramPtr getChromatogramById(int id);
  };
}

#endif
