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

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <string>
#include <vector>

namespace OpenMS
{
  const double BinaryDataEncoder_default_numpressErrorTolerance = .0001; // 1/100th of one percent
   
  /**
    @brief Class to encode and decode data encoded with MSNumpress

    MSNumpress supports three encoding schemata: 
      - Linear (MS:1002312, MS-Numpress linear prediction compression)
      - Pic (MS:1002313, MS-Numpress positive integer compression)
      - Slof (MS:1002314, MS-Numpress short logged float compression)

  */
  class OPENMS_DLLAPI MSNumpressCoder
  {

public:

    enum NumpressCompression { NONE, LINEAR, PIC, SLOF, SIZE_OF_NUMPRESSCOMPRESSION };
    /// Names of compression schemes
    static const std::string NamesOfNumpressCompression[SIZE_OF_NUMPRESSCOMPRESSION];

    /**
      @brief Configuration class for MSNumpress

      Contains configuration options for ms numpress
    */
    struct OPENMS_DLLAPI NumpressConfig 
    {
      double numpressFixedPoint;  /// fixed point for numpress algorithms
      double numpressErrorTolerance;  /// check error tolerance after encoding, guarantee abs(1.0-(encoded/decoded)) <= this, 0=do not guarantee anything
      NumpressCompression np_compression; /// which compression schema to use
      bool estimate_fixed_point; /// whether to estimate the fixed point or use the one proved with numpressFixedPoint
      double linear_fp_mass_acc; /// desired mass accuracy for linear encoding (-1 no effect, use 0.0001 for 0.2 ppm accuracy @ 500 m/z)

      NumpressConfig () :
        numpressFixedPoint(0.0),
        numpressErrorTolerance(BinaryDataEncoder_default_numpressErrorTolerance),
        np_compression(NONE),
        estimate_fixed_point(false),
        linear_fp_mass_acc(-1)
      {
      }

      /**
        @brief set compression using a string mapping to enum NumpressCompression.

        @param compression A string from NamesOfNumpressCompression[]

        @throws Exception::InvalidParameter if compression is unknown.
      */
      void setCompression(const std::string& compression)
      {
        const std::string* match = std::find(NamesOfNumpressCompression, NamesOfNumpressCompression + SIZE_OF_NUMPRESSCOMPRESSION, compression);
        
        if (match == NamesOfNumpressCompression + SIZE_OF_NUMPRESSCOMPRESSION) // == end()
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value '" + compression + "' is not a valid Numpress compression scheme.");
        }
        
        np_compression = (NumpressCompression)std::distance(NamesOfNumpressCompression, match);
      }

    };

    /// default constructor
    MSNumpressCoder() {}

    /// Destructor
    virtual ~MSNumpressCoder() {}

    /**
     * @brief Encodes a vector of floating point numbers into a Base64 string using numpress
     *
     * This code is obtained from the proteowizard implementation
     * ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest).
     *
     * This function will first apply the numpress encoding to the data, then
     * encode the result in base64 (with optional zlib compression before
     * base64 encoding).
     *
     * @note In case of error, result string is empty
     *
     * @param in The vector of floating point numbers to be encoded
     * @param result The resulting string
     * @param zlib_compression Whether to apply zlib compression after numpress compression
     * @param config The numpress configuration defining the compression strategy
     *
    */
    void encodeNP(const std::vector<double> & in, String & result,
        bool zlib_compression, const NumpressConfig & config);

    /// encodeNP from a float (convert first to double)
    void encodeNP(const std::vector<float> & in, String & result,
        bool zlib_compression, const NumpressConfig & config);

    /**
     * @brief Decodes a Base64 string to a vector of floating point numbers using numpress
     *
     * This code is obtained from the proteowizard implementation
     * ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest).
     *
     * This function will first decode the input base64 string (with optional
     * zlib decompression after decoding) and then apply numpress decoding to
     * the data.
     *
     * @param in The base64 encoded string
     * @param out The resulting vector of doubles
     * @param zlib_compression Whether to apply zlib de-compression before numpress de-compression
     * @param config The numpress configuration defining the compression strategy
     *
     * @throw throws Exception::ConversionError if the string cannot be converted
     *
    */
    void decodeNP(const String & in, std::vector<double> & out,
        bool zlib_compression, const NumpressConfig & config);

    /**
     * @brief Encode the vector in to the result string using numpress (unsafe)
     *
     * @note In case of error, result is given back unmodified
     *
     * This performs the raw numpress encoding on a set of data and does no
     * Base64 encoding on the result. Therefore the result string is likely
     * *unsafe* to handle and is basically a byte container.
     *
     * Please use the safe versions above unless you need access to the raw
     * byte arrays.
     *
     * @param in The vector of floating point numbers to be encoded
     * @param result The resulting string
     * @param config The numpress configuration defining the compression strategy
     *
    */
    void encodeNPRaw(const std::vector<double> & in, String & result, const NumpressConfig & config);

    /**
     * @brief Decode the (not necessary null terminated) string in to the result vector out (unsafe)
     *
     * @note that the string in should *only* contain the data and _no_ extra
     * null terminating byte (unless of course the last data byte is null)
     *
     * This performs the raw numpress decoding on a raw byte array (not Base64
     * encoded). Therefore the input string is likely *unsafe* to handle and is
     * basically a byte container.
     *
     * Please use the safe versions above unless you only have the raw byte
     * arrays.
     *
     * @param in The base64 encoded string
     * @param out The resulting vector of doubles
     * @param config The numpress configuration defining the compression strategy
     *
     * @throw throws Exception::ConversionError if the string cannot be converted
     *
    */
    void decodeNPRaw(const std::string & in, std::vector<double> & out, const NumpressConfig & config);

private:

    void decodeNPInternal_(const unsigned char* in, size_t in_size, std::vector<double>& out, const NumpressConfig & config);
  };

} //namespace OpenMS


