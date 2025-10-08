// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h> // StringList
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

#include <map>

namespace OpenMS
{
  /**
    @brief File adapter for MzML files

    This implementation does currently not support the whole functionality of MzML.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    MzMLFile();
    ///Destructor
    ~MzMLFile() override;

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions&);

    /**
      @brief Loads a map from a MzML file. Spectra and chromatograms are sorted by default (this can be disabled using PeakFileOptions).

      @param filename The filename with the data
      @param map Is an MSExperiment

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, PeakMap& map);

    /**
      @brief Loads a map from a MzML file stored in a buffer (in memory).

      @param[in] buffer The buffer with the data (i.e. string with content of an mzML file)
      @param[out] map Is an MSExperiment

      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadBuffer(const std::string& buffer, PeakMap& map);

    /**
      @brief Only count the number of spectra and chromatograms from a file

      This method honors PeakOptions (if specified) for spectra, i.e. only spectra within the specified
      RT range and MS levels are counted.
      If PeakOptions have no filters set (the default), then spectra and chromatogram counts
      are taken from the counts attribute of the spectrumList/chromatogramList tags (the
      parsing skips all intermediate data and ends as soon as both counts are available).

    */
    void loadSize(const String & filename, Size& scount, Size& ccount);

    /**
      @brief Stores a map in an MzML file.

      @p map has to be an MSExperiment or have the same interface.

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const PeakMap& map) const;

    /**
      @brief Stores a map in an output string.

      @p output An empty string to store the result
      @p map has to be an MSExperiment
    */
    void storeBuffer(std::string & output, const PeakMap& map) const;

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer.

      The result will not be stored directly but is available through the
      events triggered by the parser and caught by the provided IMSDataConsumer
      object.

      This function should be used if processing and storage of the result can
      be performed directly in the provided IMSDataConsumer object.

      @note Transformation can be speed up by setting skip_full_count which
      does not require a full first pass through the file to compute the
      correct number of spectra and chromatograms in the input file.

      @param filename_in Filename of input mzML file to transform
      @param consumer Consumer class to operate on the input filename (implementing a transformation)
      @param skip_full_count Whether to skip computing the correct number of spectra and chromatograms in the input file
      @param skip_first_pass Skip first file parsing pass, which hands only meta-data (number of spectra/chroms and experimental settings) to the consumer
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false, bool skip_first_pass = false);

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer

      The result will be stored in the provided map.

      This function should be used if a specific pre-processing should be
      applied to the data before storing them in a map (e.g. if data-reduction
      should be applied to the data before loading all data into memory).

      @param filename_in Filename of input mzML file to transform
      @param consumer Consumer class to operate on the input filename (implementing a transformation)
      @param map Map to store the resulting spectra and chromatograms
      @param skip_full_count Whether to skip computing the correct number of spectra and chromatograms in the input file
      @param skip_first_pass Skip first file parsing pass, which hands only meta-data (number of spectra/chroms and experimental settings) to the consumer
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, PeakMap& map, bool skip_full_count = false, bool skip_first_pass = false);

    /**
      @brief Checks if a file validates against the XML schema.

      @exception Exception::FileNotFound is thrown if the file cannot be found.
    */
    bool isValid(const String& filename, std::ostream& os = std::cerr);

    /**
      @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

      @param filename File name of the file to be checked.
      @param errors Errors during the validation are returned in this output parameter.
      @param warnings Warnings during the validation are returned in this output parameter.

      @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings);

    /**
      @brief Checks if a file is an indexed MzML file or not.

      @param filename File name of the file to be checked.

      @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool hasIndex(const String& filename);

    struct SpecInfo
    {
      Size count_centroided = 0;
      Size count_profile = 0;
      Size count_unknown = 0;
    };
      
    /**
       @brief Check type of spectra based on their metadata (if available) or by inspecting the peaks itself.
       
       By default, only the first @p first_n_spectra_only, which are NOT 'unknown' are checked to save time.
       The current PeakFileOptions, e.g. which MS-level to read/skip, are honored, e.g. skipped spectra do not count
       towards @p first_n_spectra_only.
       
       You can use this function to estimate the spectrum type, but it should be done for each MS-level separately.
       Otherwise you might get mixed (PROFILE+CENTROIDED) results.
       
       @param filename File name of the mzML file to be checked
       @param first_n_spectra_only Only inspect this many spectra (UNKNOWN spectra do not count) and then end parsing the file
       
       @return Map of MS level to counts (centroided, profile, unknown)
       
       @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    std::map<UInt, SpecInfo> getCentroidInfo(const String& filename, const Size first_n_spectra_only = 10);

protected:

    /// Perform first pass through the file and retrieve the meta-data to initialize the consumer
    void transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count);

    /// Safe parse that catches exceptions and handles them accordingly
    void safeParse_(const String & filename, Internal::XMLHandler * handler);

private:

    /// Options for loading / storing
    PeakFileOptions options_;

    /// Location of indexed mzML schema
    String indexed_schema_location_;
  };

} // namespace OpenMS
