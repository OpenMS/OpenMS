// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>

#include <string>
#include <xercesc/dom/DOMNode.hpp>

#include <OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

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

    bool skip_xml_checks_; ///< Whether to skip some XML checks (e.g. removing whitespace inside base64 arrays) and be fast instead
      
    typedef Internal::MzMLHandlerHelper::BinaryData BinaryData;

    /**
      @brief decode binary data

      @todo Duplicated code from MzMLHandler, need to clean up see MzMLHandler::fillData_() 

    */
    OpenMS::Interfaces::SpectrumPtr decodeBinaryDataSpectrum_(std::vector<BinaryData> & data);

    void decodeBinaryDataMSSpectrum_(std::vector<BinaryData>& data, OpenMS::MSSpectrum& s);

    void decodeBinaryDataMSChrom_(std::vector<BinaryData>& data, OpenMS::MSChromatogram& c);

    /**
      @brief decode binary data

      @todo Duplicated code from MzMLHandler, need to clean up see MzMLHandler::fillData_() 

    */
    OpenMS::Interfaces::ChromatogramPtr decodeBinaryDataChrom_(std::vector<BinaryData> & data);

    /**
      @brief Convert a single DOMNode of type binaryDataArray to BinaryData object.

      This function will extract the data from a xerces DOMNode which points to
      a binaryDataArray tag and store the result as a BinaryData object. The
      result will be appended to the data vector.

      @param in DOMNode of type binaryDataArray
      @param data Binary data extracted from the string
    */
    void handleBinaryDataArray_(xercesc::DOMNode * indexListNode, std::vector<BinaryData>& data);

    /**
      @brief Extract data from a string containing multiple <binaryDataArray> tags.

      This may be a string from <spectrum> to </spectrum> or <chromatogram> to
      </chromatogram> tag which contains one or more <binaryDataArray>. These
      XML tags need to conform to the mzML standard. The function will return a
      vector with all binary data found in the string in the binaryDataArray
      tags.

      @param in Input string containing the raw XML
      @param data Binary data extracted from the string

      @pre in must have <spectrum> or <chromatogram> as root element.

    */
    std::string domParseString_(const std::string& in, std::vector<BinaryData>& data);

  public:

    explicit MzMLSpectrumDecoder(bool skip_xml_checks = false) :
      skip_xml_checks_(skip_xml_checks)
    {}

    /**
      @brief Extract data from a string which contains a full mzML spectrum.

      Extracts data from the input string which is expected to contain exactly
      one <spectrum> tag (from <spectrum> to </spectrum>). This function will
      extract the contained binaryDataArray and provide the result as Spectrum.

      @param in Input string containing the raw XML
      @param sptr Resulting spectrum

      @pre in must have <spectrum> as root element.

    */
    void domParseSpectrum(const std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr);

    /**
      @brief Extract data from a string which contains a full mzML spectrum.

      Extracts data from the input string which is expected to contain exactly
      one <spectrum> tag (from <spectrum> to </spectrum>). This function will
      extract the contained binaryDataArray and provide the result as Spectrum.

      @param in Input string containing the raw XML
      @param s Resulting spectrum

      @pre in must have <spectrum> as root element.

    */
    void domParseSpectrum(const std::string& in, MSSpectrum& s);

    /**
      @brief Extract data from a string which contains a full mzML chromatogram.

      Extracts data from the input string which is expected to contain exactly
      one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
      function will extract the contained binaryDataArray and provide the
      result as Chromatogram.

      @param in Input string containing the raw XML
      @param c Resulting chromatogram

      @pre in must have <chromatogram> as root element.
    */
    void domParseChromatogram(const std::string& in, MSChromatogram& c);

    /**
      @brief Extract data from a string which contains a full mzML chromatogram.

      Extracts data from the input string which is expected to contain exactly
      one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
      function will extract the contained binaryDataArray and provide the
      result as Chromatogram.

      @param in Input string containing the raw XML
      @param cptr Resulting chromatogram

      @pre in must have <chromatogram> as root element.
    */
    void domParseChromatogram(const std::string& in, OpenMS::Interfaces::ChromatogramPtr & cptr);

    /// Whether to skip some XML checks (e.g. removing whitespace inside base64 arrays) and be fast instead
    void setSkipXMLChecks(bool only);
  };
}


