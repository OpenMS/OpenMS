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

#pragma once

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <vector>

namespace OpenMS
{
  namespace Internal
  {

    /**@brief Helper for mzML file format
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

      /**@brief Representation for binary data in mzML
       *
       * Represents data in the <binaryDataArray> tag
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
                                        String indent = "",
                                        bool use_numpress = false);

      /**
        @brief Write the indexed mzML footer the appropriate compression term given the PeakFileOptions and the NumpressConfig

        @param os The output stream
        @param options The PeakFileOptions used for writing
        @param spectra_offsets Binary offsets of <spectrum> tags
        @param chromatograms_offsets Binary offsets of <chromatogram> tags

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


