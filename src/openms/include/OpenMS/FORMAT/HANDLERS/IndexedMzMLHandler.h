// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/INTERFACES/ISpectrumAccess.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <string>
#include <fstream>
#include <unordered_map>

namespace OpenMS
{

namespace Internal
{

  /**
    @brief A low-level class to read an indexedmzML file.

    This class provides low-level access to the underlying data structures, if
    you simply want to read an indexed mzML file you probably want to use
    IndexedMzMLFileLoader instead.

    This class implements access to an indexedmzML file and the contained spectra
    and chromatogram data through the getSpectrumById and getChromatogramById
    functions. It thus allows random access to spectra and chromatograms data
    without having to read the whole file into memory. It does not provide the
    same interface as MSExperiment, if this is desired, please use
    IndexedMzMLFileLoader and OnDiscMSExperiment.

    Internally, it uses the IndexedMzMLDecoder for initial parsing and
    extracting all the offsets of the &lt;chromatogram&gt; and &lt;spectrum&gt; tags. These
    offsets are stored as members of this class as well as the offset to the &lt;indexList&gt; element

    @note This implementation is @a not thread-safe since it keeps internally a
    single file access pointer which it moves when accessing a specific
    data item. The caller is responsible to ensure that access is performed
    atomically.

  */
  class OPENMS_DLLAPI IndexedMzMLHandler
  {
    /// Name of the file
    String filename_;
    /// Binary offsets to all spectra
    std::vector< std::streampos > spectra_offsets_;
    /// Mapping of spectra native ids to offsets
    std::unordered_map< std::string, Size > spectra_native_ids_;
    /// Binary offsets to all chromatograms
    std::vector< std::streampos > chromatograms_offsets_;
    /// Mapping of chromatogram native ids to offsets
    std::unordered_map< std::string, Size > chromatograms_native_ids_;
    /// offset to the \<indexList\> element
    std::streampos index_offset_;
    /// Whether spectra are written before chromatograms in this file
    bool spectra_before_chroms_;
    /// The current filestream (opened by openFile)
    std::ifstream filestream_;
    /// Whether parsing the indexedmzML file was successful
    bool parsing_success_;
    /// Whether to skip XML checks
    bool skip_xml_checks_;

    /**
      @brief Try to parse the footer of the indexedmzML

      Upon success, the chromatogram and spectra offsets will be populated and
      parsing_success_ will be set to true.

      @note You *need* to check getParsingSuccess after calling this!
    */
    void parseFooter_();

    std::string getChromatogramById_helper_(int id);

    std::string getSpectrumById_helper_(int id);

    public:

    /**
      @brief Default constructor
    */
    IndexedMzMLHandler();

    /**
      @brief Constructor

      Tries to parse the file, success can be checked with getParsingSuccess()
    */
    explicit IndexedMzMLHandler(const String& filename);

    /// Copy constructor
    IndexedMzMLHandler(const IndexedMzMLHandler& source);

    /// Destructor
    ~IndexedMzMLHandler();

    /**
      @brief Open a file

      Tries to parse the file, success can be checked with getParsingSuccess()
    */
    void openFile(const String& filename);

    /**
      @brief Returns whether parsing was successful

      @note Callable after openFile or the constructor using a filename
      @note It is invalid to call getSpectrumById or getChromatogramById if this function returns false

      @return Whether the parsing of the file was successful (if false, the
      file most likely was not an indexed mzML file)
    */
    bool getParsingSuccess() const;

    /// Returns the number of spectra available
    size_t getNrSpectra() const;

    /// Returns the number of chromatograms available
    size_t getNrChromatograms() const;

    /**
      @brief Retrieve the raw data for the spectrum at position "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id is not within [0, getNrSpectra()-1]

      @return The spectrum at position id
    */
    OpenMS::Interfaces::SpectrumPtr getSpectrumById(int id);

    /**
      @brief Retrieve the raw data for the spectrum at position "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id is not within [0, getNrSpectra()-1]

      @return The spectrum at position id
    */
    const OpenMS::MSSpectrum getMSSpectrumById(int id);

    /**
      @brief Retrieve the raw data for the spectrum with native id "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id cannot be found

      @param id The spectrum native id
      @param s The spectrum to be used and filled with data
    */
    void getMSSpectrumByNativeId(const std::string& id, OpenMS::MSSpectrum& s);

    /**
      @brief Retrieve the raw data for the spectrum at position "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id is not within [0, getNrSpectra()-1]

      @param id The spectrum id
      @param s The spectrum to be used and filled with data
    */
    void getMSSpectrumById(int id, OpenMS::MSSpectrum& s);

    /**
      @brief Retrieve the raw data for the chromatogram at position "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id is not within [0, getNrChromatograms()-1]

      @return The chromatogram at position id
    */
    OpenMS::Interfaces::ChromatogramPtr getChromatogramById(int id);

    /**
      @brief Retrieve the raw data for the chromatogram at position "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id is not within [0, getNrChromatograms()-1]

      @return The chromatogram at position id
    */
    const OpenMS::MSChromatogram getMSChromatogramById(int id);

    /**
      @brief Retrieve the raw data for the chromatogram with native id "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id cannot be found

      @param id The chromatogram native id
      @param c The chromatogram to be used and filled with data
    */
    void getMSChromatogramByNativeId(const std::string& id, OpenMS::MSChromatogram& c);

    /**
      @brief Retrieve the raw data for the chromatogram at position "id"

      @throw Exception if getParsingSuccess() returns false
      @throw Exception if id is not within [0, getNrChromatograms()-1]

      @param id The chromatogram id
      @param c The chromatogram to be used and filled with data
    */
    void getMSChromatogramById(int id, OpenMS::MSChromatogram& c);

    /// Whether to skip some XML checks (removing whitespace from base64 arrays) and be fast instead
    void setSkipXMLChecks(bool skip)
    {
      skip_xml_checks_ = skip;
    }

  };
}
}

