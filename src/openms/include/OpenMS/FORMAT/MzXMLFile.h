// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{
  class String;

  /**
      @brief File adapter for MzXML 3.1 files

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
	typedef PeakMap MapType;

public:
    ///Default constructor
    MzXMLFile();
    ///Destructor
    ~MzXMLFile() override;

    /// Mutable access to the options for loading/storing
    PeakFileOptions & getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions & getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions &);

    /**
        @brief Loads a map from a MzXML file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, MapType & map);

    /**
        @brief Stores a map in a MzXML file.

        @p map has to be a MSExperiment or have the same interface.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const MapType & map) const;

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
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false);

    /**
      @brief Transforms a map while loading using the supplied MSDataConsumer

      The result will be stored in the provided map.

      This function should be used if a specific pre-processing should be
      applied to the data before storing them in a map (e.g. if data-reduction
      should be applied to the data before loading all data into memory).

      @param filename_in Filename of input mzXML file to transform 
      @param consumer Consumer class to operate on the input filename (implementing a transformation)
      @param map Map to store the resulting spectra and chromatograms
      @param skip_full_count Whether to skip computing the correct number of spectra and chromatograms in the input file 
    */
    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, MapType& map, bool skip_full_count = false);

protected:

    /// Perform first pass through the file and retrieve the meta-data to initialize the consumer
    void transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count);

private:

    PeakFileOptions options_;
  };
} // namespace OpenMS

