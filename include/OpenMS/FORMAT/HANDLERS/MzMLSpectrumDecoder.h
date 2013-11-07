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

#ifndef OPENMS_FORMAT_HANDLERS_MZMLSPECTRUMDECODER_H
#define OPENMS_FORMAT_HANDLERS_MZMLSPECTRUMDECODER_H

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <string>
#include <xercesc/dom/DOMNode.hpp>

namespace OpenMS
{

  /**
    @brief A class to decode input strings that contain an mzML chromatogram or
    spectrum tag.

    It uses xercesc to parse a string containing either a exactly one mzML
    spectrum or chromatogram (from <chromatogram> to </chromatogram> or
    <spectrum> to </spectrum> tag). It returns the data contained in the
    binaryDataArray for Intensity / mass-to-charge or Intensity / time.

  */
  class OPENMS_DLLAPI MzMLSpectrumDecoder
  {
  protected:
    /// Spectrum representation
    struct BinaryData
    {
      String base64;
      enum {PRE_NONE, PRE_32, PRE_64} precision;
      Size size;
      bool compression;
      enum {DT_NONE, DT_FLOAT, DT_INT, DT_STRING} data_type;
      std::vector<Real> floats_32;
      std::vector<DoubleReal> floats_64;
      std::vector<Int32> ints_32;
      std::vector<Int64> ints_64;
      std::vector<String> decoded_char;
      MetaInfoDescription meta;
    };

    /**
      @brief Decode base64 arrays

      @TODO Duplicated code from MzMLHandler, need to clean up
      see void MzMLHandler<MapType>::fillData_() first 60 LOC

    */
    void decode64arrays(std::vector<BinaryData> & data_);

    /**
      @brief compute data properties

      @TODO Duplicated code from MzMLHandler, need to clean up
      see void MzMLHandler<MapType>::fillData_() 

    */
    void computeDataProperties_(std::vector<BinaryData> & data_, 
      bool& precision_64, SignedSize& index, String index_name);
    
    /**
      @brief decode binary data

      @TODO Duplicated code from MzMLHandler, need to clean up
      see void MzMLHandler<MapType>::fillData_() 

    */
    OpenMS::Interfaces::SpectrumPtr decodeBinaryData(std::vector<BinaryData> & data_);

    /**
      @brief decode binary data

      @TODO Duplicated code from MzMLHandler, need to clean up
      see void MzMLHandler<MapType>::fillData_() 

    */
    OpenMS::Interfaces::ChromatogramPtr decodeBinaryDataChrom(std::vector<BinaryData> & data_);

    /**
      @brief Handle (parent_tag == "binaryDataArray") cv term

      @TODO Duplicated code from MzMLHandler, need to clean up
      see void MzMLHandler<MapType>::handleCVParam_(...) 

    */
    void handleCVParam(std::vector<BinaryData> & data_, 
        const String& accession, const String& value, const String& name);

    /**
      @brief Convert a single DOMNode of type binaryDataArray to BinaryData object.

      This function will extract the data from a xerces DOMNode which points to
      a binaryDataArray tag and store the result as a BinaryData object. The
      result will be appended to the data_ vector.

      @param in DOMNode of type binaryDataArray
      @param data_ Binary data extracted from the string
    */
    void handleBinaryDataArray(xercesc::DOMNode * indexListNode, std::vector<BinaryData>& data_);

    /**
      @brief Extract data from a string containing multiple binaryDataArray tags.

      This may be a string from <spectrum> to </spectrum> or <chromatogram> to
      </chromatogram> tag which contains one or more binaryDataArray.  The
      function will return a vector with all binary data found in the string in
      the binaryDataArray tags.

      @param in Input string containing the raw XML
      @param data_ Binary data extracted from the string
    */
    void domParseString(const std::string& in, std::vector<BinaryData>& data_);

  public:

    /**
      @brief Extract data from a string which contains a full mzML spectrum.

      Extracts data from the input string which is expected to contain exactly
      one <spectrum> tag (from <spectrum> to </spectrum>). This function will
      extract the contained binaryDataArray and provide the result as Spectrum.
    */
    void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr);

    /**
      @brief Extract data from a string which contains a full mzML chromatogram.

      Extracts data from the input string which is expected to contain exactly
      one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
      function will extract the contained binaryDataArray and provide the
      result as Chromatogram.
    */
    void domParseChromatogram(std::string& in, OpenMS::Interfaces::ChromatogramPtr & sptr);

  };
}

#endif
