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

#ifndef OPENMS_FORMAT_NUMPRESSCODER_H
#define OPENMS_FORMAT_NUMPRESSCODER_H

#include <OpenMS/FORMAT/Base64.h>

namespace OpenMS
{
  const double BinaryDataEncoder_default_numpressErrorTolerance = .0001; // 1/100th of one percent
   
  /**
    @brief Class to encode and decode data encoded with MSNumpress

    MSNumpress supports three encoding schemas: 
      - Linear (MS:1002312, MS-Numpress linear prediction compression)
      - Pic (MS:1002313, MS-Numpress positive integer compression)
      - Slof (MS:1002314, MS-Numpress short logged float compression)

  */
  class OPENMS_DLLAPI MSNumpressCoder
  {

public:

    enum NumpressCompression { NONE, LINEAR, PIC, SLOF };

    /**
      @brief Configuration class for MSNumpress

      Contains configuration options for ms numpress
    */
    struct NumpressConfig 
    {
      double numpressFixedPoint;  /// fixed point for numpress algorithms
      double numpressErrorTolerance;  /// check error tolerance after encoding, guarantee abs(1.0-(encoded/decoded)) <= this, 0=do not guarantee anything
      NumpressCompression np_compression; /// which compression schema to use
      bool estimate_fixed_point; /// whether to estimate the fixed point or use the one proved with numpressFixedPoint

      NumpressConfig () :
        numpressFixedPoint(0.0),
        numpressErrorTolerance(BinaryDataEncoder_default_numpressErrorTolerance),
        np_compression(NONE),
        estimate_fixed_point(false)
      {
      }

    };

    /// default constructor
    MSNumpressCoder() {}

    /// Destructor
    virtual ~MSNumpressCoder() {}

    /**
        @brief Encodes a vector of floating point numbers into a Base64 string using numpress

        This code is obtained from the proteowizard implementation
        ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp

        adapted by Hannes Roest

        In case of error, result is empty
    */
    void encodeNP(std::vector<double> & in, String & result,
        bool zlib_compression, NumpressConfig config)
    {
      result.clear();
      encodeNP_(in, result, config);
      if (result.empty())
      {
        return;
      }

      // Encode in base64 and compress
      std::vector<String> tmp;
      tmp.push_back(result);
      base64coder.encodeStrings(tmp, result, zlib_compression, false);
    }

    /// encodeNP from a float (convert first to double)
    void encodeNP(std::vector<float> & in, String & result,
        bool zlib_compression, NumpressConfig config)
    {
      std::vector<float> dvector(in);
      encodeNP(in, result, zlib_compression, config);
    }

    /**
        @brief Decodes a Base64 string to a vector of floating point numbers using numpress

        This code is obtained from the proteowizard implementation
        ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp

        adapted by Hannes Roest
    */
    void decodeNP(const String & in, std::vector<double> & out,
        bool zlib_compression, NumpressConfig config)
    {
      QByteArray base64_uncompressed;
      base64coder.decodeSingleString(in, base64_uncompressed, zlib_compression);

      // Create a temporary string (*not* null-terminated) to hold the data
      std::string tmpstring(base64_uncompressed.constData(), base64_uncompressed.size());
      decodeNP_(tmpstring, out, config);

      // NOTE: it is possible (and likely faster) to call directly the const
      // unsigned char * function but this would make it necessary to deal with
      // reinterpret_cast ugliness here ... 
      //
      // decodeNP_internal_(reinterpret_cast<const unsigned char*>(base64_uncompressed.constData()), base64_uncompressed.size(), out, config);
    }

private:

    /**
        @brief Encode the vector in to the result string

        In case of error, result is given back unmodified
    */
    void encodeNP_(std::vector<double> & in, String & result, NumpressConfig config);

    /**
        @brief Decode the (not necessary null terminated) string in to the result vector out

        In case of error, an exception is thrown

        @note that the string in should *only* contain the data and _no_ extra
        null terminating byte (unless of course the last data byte is null)
    */
    void decodeNP_(std::string in, std::vector<double> & out, NumpressConfig config);

    void decodeNP_internal_(const unsigned char* in, size_t in_size, std::vector<double> & out, NumpressConfig config);

    Base64 base64coder;
  };

} //namespace OpenMS

#endif /* OPENMS_FORMAT_MSNUMPRESSCODER_H */
