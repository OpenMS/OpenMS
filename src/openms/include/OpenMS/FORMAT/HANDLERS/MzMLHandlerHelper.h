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

#ifndef OPENMS_FORMAT_HANDLERS_MZMLHANDLERHELPER_H
#define OPENMS_FORMAT_HANDLERS_MZMLHANDLERHELPER_H

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
        String base64;
        enum {PRE_NONE, PRE_32, PRE_64} precision;
        Size size;
        bool compression; // zlib compression
        enum {DT_NONE, DT_FLOAT, DT_INT, DT_STRING} data_type;
        std::vector<Real> floats_32;
        std::vector<DoubleReal> floats_64;
        std::vector<Int32> ints_32;
        std::vector<Int64> ints_64;
        std::vector<String> decoded_char;
        MetaInfoDescription meta;
        MSNumpressCoder::NumpressCompression np_compression;
      };

      /**
        @brief Returns the appropriate compression term given the PeakFileOptions and the NumpressConfig  
      */
      static String getCompressionTerm_(const PeakFileOptions& opt, MSNumpressCoder::NumpressConfig np_compression, bool use_numpress = false);

      /**
        @brief Write the mzML footer the appropriate compression term given the PeakFileOptions and the NumpressConfig  
      */
      static void writeFooter_(std::ostream& os, const PeakFileOptions& options_,
        std::vector< std::pair<std::string, long> > & spectra_offsets,
        std::vector< std::pair<std::string, long> > & chromatograms_offsets
      );

      static void decodeBase64Arrays(std::vector<BinaryData> & data_);

      static void computeDataProperties_(std::vector<BinaryData>& data_, bool& precision_64, SignedSize& index, String index_name);

      static bool handleBinaryDataArrayCVParam(std::vector<BinaryData>& data_,
        const String& accession, const String& value, const String& name);
    };


  } // namespace Internal
} // namespace OpenMS

#endif
