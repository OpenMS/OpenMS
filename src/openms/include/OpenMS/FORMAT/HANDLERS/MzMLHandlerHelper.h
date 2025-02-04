// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <vector>

namespace OpenMS
{
  namespace Internal
  {

    /**
     * @brief Helper for mzML file format
     *
     * This class provides common structures and re-useable helper functions
     * for parsing the mzML format. These are mainly used by MzMLHandler and MzMLSpectrumDecoder.
     *
     **/
    class OPENMS_DLLAPI MzMLHandlerHelper
    {

      /// Also display some warning message when appropriate (see XMLHandler)
      static void warning(int mode, const String & msg, UInt line = 0, UInt column = 0);

    public:

      /**
       * @brief Representation for binary data in mzML
       *
       * Represents data in the `<binaryDataArray>` tag
       *
       **/
      struct BinaryData
      {
        // ordered by size (alignment) and cache hotness in 'decode'

        enum {
          PRE_NONE, ///< unknown precision
          PRE_32,   ///< 32bit precision
          PRE_64    ///< 64bit precision
        } precision;

        enum {
          DT_NONE,    ///< unknown data type
          DT_FLOAT,   ///< float data type
          DT_INT,     ///< integer data type
          DT_STRING   ///< string data type
        } data_type;

        MSNumpressCoder::NumpressCompression np_compression; ///< numpress options

        bool compression; ///< zlib compression
        double unit_multiplier; ///< multiplier for unit (e.g. 60 for minutes)

        String base64; ///< Raw data in base64 encoding
        Size size; ///< Raw data length
        std::vector<float> floats_32;
        std::vector<double> floats_64;
        std::vector<Int32> ints_32;
        std::vector<Int64> ints_64;
        std::vector<String> decoded_char;

        MetaInfoDescription meta; ///< Meta data description

        /// Constructor
        BinaryData() :
          precision(PRE_NONE),
          data_type(DT_NONE),
          np_compression(),
          compression(false),
          unit_multiplier(1.0),
          base64(),
          size(0),
          floats_32(),
          floats_64(),
          ints_32(),
          ints_64(),
          decoded_char(),
          meta()
        {
        }

        BinaryData(const BinaryData&) = default;               // Copy constructor
        BinaryData(BinaryData&&) = default;                    // Move constructor
        BinaryData& operator=(const BinaryData&) & = default;  // Copy assignment operator
        BinaryData& operator=(BinaryData&&) & = default;       // Move assignment operator
        ~BinaryData() = default;                               // Destructor

      };

      /**
        @brief Returns the appropriate compression term given the PeakFileOptions and the NumpressConfig
      */
      static String getCompressionTerm_(const PeakFileOptions& opt,
                                        MSNumpressCoder::NumpressConfig np_compression,
                                        const String& indent = "",
                                        bool use_numpress = false);

      /**
        @brief Write the indexed mzML footer the appropriate compression term given the PeakFileOptions and the NumpressConfig

        @param os The output stream
        @param options The PeakFileOptions used for writing
        @param spectra_offsets Binary offsets of &lt;spectrum&gt; tags
        @param chromatograms_offsets Binary offsets of &lt;chromatogram&gt; tags

      */
      static void writeFooter_(std::ostream& os,
                               const PeakFileOptions& options,
                               const std::vector< std::pair<std::string, Int64> > & spectra_offsets,
                               const std::vector< std::pair<std::string, Int64> > & chromatograms_offsets);

      /**
        @brief Decode Base64 arrays and write into data_ array

        @param data_ The input and output
        @param skipXMLCheck whether to skip cleaning the Base64 arrays and remove whitespaces
      */
      static void decodeBase64Arrays(std::vector<BinaryData> & data_, const bool skipXMLCheck = false);

      /**
        @brief Identify a data array from a list.

        Given a specific array name, find it in the provided list and return its index and precision.

        @param data_ The list of data arrays
        @param precision_64 Whether the identified array has 64 bit precision
        @param index The index of the identified array
        @param index_name The name of the array to be identified
      */
      static void computeDataProperties_(const std::vector<BinaryData>& data_, bool& precision_64, SignedSize& index, const String& index_name);

      /**
        @brief Handle a given CV parameter found in a binaryDataArray tag

        Given a CV parameter, properly set the members of the last entry of
        data_, this will properly handle all terms describing precision,
        compression, name of the data and units.

        @param data_ The list of data arrays, whose last entry will be changed
        @param accession The CV accession
        @param value The CV value
        @param name The CV name
        @param unit_accession The CV unit accession (if a unit tag is present)
      */
      static bool handleBinaryDataArrayCVParam(std::vector<BinaryData>& data_,
                                               const String& accession,
                                               const String& value,
                                               const String& name,
                                               const String& unit_accession);
    };


  } // namespace Internal
} // namespace OpenMS


