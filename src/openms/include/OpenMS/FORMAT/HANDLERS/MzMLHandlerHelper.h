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

#pragma once

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <vector>

namespace OpenMS
{
  namespace Internal
  {

    class OPENMS_DLLAPI MzMLHandlerHelper
    {

      /// Also display some warning message when appropriate (see XMLHandler)
      static void warning(int mode, const String & msg, UInt line = 0, UInt column = 0);

    public:

      /// Binary data representation
      struct BinaryData
      {
        // ordered by size (alignment) and cache hotness in 'decode'
        enum {PRE_NONE, PRE_32, PRE_64} precision;
        enum { DT_NONE, DT_FLOAT, DT_INT, DT_STRING } data_type;
        MSNumpressCoder::NumpressCompression np_compression;
        bool compression; // zlib compression
        double unit_multiplier;
        String base64;
        Size size;
        std::vector<float> floats_32;
        std::vector<double> floats_64;
        std::vector<Int32> ints_32;
        std::vector<Int64> ints_64;
        std::vector<String> decoded_char;
        MetaInfoDescription meta;

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

      };

      /**
        @brief Returns the appropriate compression term given the PeakFileOptions and the NumpressConfig  
      */
      static String getCompressionTerm_(const PeakFileOptions& opt, MSNumpressCoder::NumpressConfig np_compression, String indent = "", bool use_numpress = false);

      /**
        @brief Write the mzML footer the appropriate compression term given the PeakFileOptions and the NumpressConfig  
      */
      static void writeFooter_(std::ostream& os, const PeakFileOptions& options_,
        std::vector< std::pair<std::string, long> > & spectra_offsets,
        std::vector< std::pair<std::string, long> > & chromatograms_offsets
      );

      /**
        @brief Decode Base64 arrays and write into data_ array
        
        @param data_ The input and output
        @param skipXMLCheck whether to skip cleaning the Base64 arrays and remove whitespaces 
      */
      static void decodeBase64Arrays(std::vector<BinaryData> & data_, const bool skipXMLCheck = false);

      static void computeDataProperties_(const std::vector<BinaryData>& data_, bool& precision_64, SignedSize& index, const String& index_name);

      static bool handleBinaryDataArrayCVParam(std::vector<BinaryData>& data_,
        const String& accession, const String& value, const String& name, const String& unit_accession);
    };


  } // namespace Internal
} // namespace OpenMS


