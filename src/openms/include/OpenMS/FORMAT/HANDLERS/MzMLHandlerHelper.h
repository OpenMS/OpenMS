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
#include <OpenMS/FORMAT/ZstdCompression.h>
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
      static void warning(int mode, const std::string & msg, UInt line = 0, UInt column = 0);

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

        enum PRECISION{
          PRE_NONE, ///< unknown precision
          PRE_32,   ///< 32bit precision
          PRE_64    ///< 64bit precision
        } precision;

        enum DATA_TYPE{
          DT_NONE,    ///< unknown data type
          DT_FLOAT,   ///< float data type
          DT_INT,     ///< integer data type
          DT_STRING   ///< string data type
        } data_type;

        MSNumpressCoder::NumpressCompression np_compression; ///< numpress options

        bool compression; ///< zlib compression
        bool zstd_compression; ///< zstd compression (mutually exclusive with zlib compression)
        ZstdCompression::ByteTransform zstd_transform; ///< transform applied before zstd compression (only used if zstd_compression is set)
        double unit_multiplier; ///< multiplier for unit (e.g. 60 for minutes)

        std::string base64; ///< Raw data in base64 encoding
        Size size; ///< Raw data length
        std::vector<float> floats_32;
        std::vector<double> floats_64;
        std::vector<Int32> ints_32;
        std::vector<Int64> ints_64;
        std::vector<std::string> decoded_char;

        MetaInfoDescription meta; ///< Meta data description

        /// Constructor
        BinaryData() :
          precision(PRE_NONE),
          data_type(DT_NONE),
          np_compression(),
          compression(false),
          zstd_compression(false),
          zstd_transform(ZstdCompression::ByteTransform::NONE),
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

        @param[in] opt The PeakFileOptions used for writing (zlib or zstd compression)
        @param[in] np_compression The numpress configuration of the array
        @param[in] indent Indentation prepended to the term
        @param[in] use_numpress Whether numpress compression is applied to the array
        @param[in] zstd_byte_shuffle Whether the array is byte-shuffled before zstd compression (only used if zstd compression is enabled and numpress is not used)
      */
      static std::string getCompressionTerm_(const PeakFileOptions& opt,
                                        MSNumpressCoder::NumpressConfig np_compression,
                                        const std::string& indent = "",
                                        bool use_numpress = false,
                                        bool zstd_byte_shuffle = true);

      /**
        @brief Encode a numeric array as Base64 string, applying the compression selected in @p opt

        If zstd compression is enabled in @p opt, the array is byte-shuffled and compressed with
        zstd (MS:1003781), otherwise it is zlib-compressed if requested (MS:1000574). The data is
        always stored in little-endian byte order.

        @param[in,out] in The data to encode (may be modified, i.e. endianized)
        @param[in] opt The PeakFileOptions used for writing
        @param[out] out The resulting Base64 string
      */
      static void encodeNumericArray(std::vector<float>& in, const PeakFileOptions& opt, std::string& out);
      /// @copydoc encodeNumericArray(std::vector<float>&, const PeakFileOptions&, std::string&)
      static void encodeNumericArray(std::vector<double>& in, const PeakFileOptions& opt, std::string& out);
      /// @copydoc encodeNumericArray(std::vector<float>&, const PeakFileOptions&, std::string&)
      static void encodeNumericArray(std::vector<Int32>& in, const PeakFileOptions& opt, std::string& out);
      /// @copydoc encodeNumericArray(std::vector<float>&, const PeakFileOptions&, std::string&)
      static void encodeNumericArray(std::vector<Int64>& in, const PeakFileOptions& opt, std::string& out);

      /**
        @brief Encode a string array (null-terminated strings) as Base64 string, applying the compression selected in @p opt

        If zstd compression is enabled in @p opt, the data is compressed with zstd (MS:1003780),
        otherwise it is zlib-compressed if requested (MS:1000574).

        @param[in] in The strings to encode
        @param[in] opt The PeakFileOptions used for writing
        @param[out] out The resulting Base64 string
      */
      static void encodeStringArray(const std::vector<std::string>& in, const PeakFileOptions& opt, std::string& out);

      /**
        @brief Numpress-encode an array and encode the result as Base64 string, applying the compression selected in @p opt

        If zstd compression is enabled in @p opt, the numpress output is compressed with zstd
        (MS:1003783 - MS:1003785), otherwise it is zlib-compressed if requested (MS:1002746 - MS:1002748).

        @param[in] in The data to encode
        @param[in] opt The PeakFileOptions used for writing
        @param[in] config The numpress configuration
        @param[out] out The resulting Base64 string (empty if numpress encoding failed)
      */
      static void encodeNumpressArray(const std::vector<double>& in, const PeakFileOptions& opt, const MSNumpressCoder::NumpressConfig& config, std::string& out);
      /// @copydoc encodeNumpressArray(const std::vector<double>&, const PeakFileOptions&, const MSNumpressCoder::NumpressConfig&, std::string&)
      static void encodeNumpressArray(const std::vector<float>& in, const PeakFileOptions& opt, const MSNumpressCoder::NumpressConfig& config, std::string& out);

      /**
        @brief Write the indexed mzML footer the appropriate compression term given the PeakFileOptions and the NumpressConfig

        @param[out] os The output stream
        @param[in] options The PeakFileOptions used for writing
        @param[in] spectra_offsets Binary offsets of &lt;spectrum&gt; tags
        @param[in] chromatograms_offsets Binary offsets of &lt;chromatogram&gt; tags

      */
      static void writeFooter_(std::ostream& os,
                               const PeakFileOptions& options,
                               const std::vector< std::pair<std::string, Int64> > & spectra_offsets,
                               const std::vector< std::pair<std::string, Int64> > & chromatograms_offsets);

      /**
        @brief Decode Base64 arrays and write into data_ array

        @param[in,out] data_ The input and output
        @param[in] skipXMLCheck whether to skip cleaning the Base64 arrays and remove whitespaces
      */
      static void decodeBase64Arrays(std::vector<BinaryData> & data_, const bool skipXMLCheck = false);

      /**
        @brief Identify a data array from a list.

        Given a specific array name, find it in the provided list and return its index and precision.

        @param[in] data_ The list of data arrays
        @param[in] precision_64 Whether the identified array has 64 bit precision
        @param[in] index The index of the identified array
        @param[in] index_name The name of the array to be identified
      */
      static void computeDataProperties_(const std::vector<BinaryData>& data_, bool& precision_64, SignedSize& index, const std::string& index_name);

      /**
        @brief Handle a given CV parameter found in a binaryDataArray tag

        Given a CV parameter, properly set the members of the last entry of
        data_, this will properly handle all terms describing precision,
        compression, name of the data and units.

        @param[in,out] data_ The list of data arrays, whose last entry will be changed
        @param[in] accession The CV accession
        @param[in] value The CV value
        @param[in] name The CV name
        @param[in] unit_accession The CV unit accession (if a unit tag is present)
      */
      static bool handleBinaryDataArrayCVParam(std::vector<BinaryData>& data_,
                                               const std::string& accession,
                                               const std::string& value,
                                               const std::string& name,
                                               const std::string& unit_accession);
    };


  } // namespace Internal
} // namespace OpenMS


