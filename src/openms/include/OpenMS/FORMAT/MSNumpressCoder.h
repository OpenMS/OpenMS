// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <algorithm>
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

    Note that the linear compression scheme only makes sense for monotonically
    increasing data (such as retention time and m/z) that is often equally
    spaced. Pic compression only makes sense for positive integers as all data
    will be rounded to the nearest integer. Slof makes sense for all other data
    (such as non-integer intensity values).

    For more information on the compression schemata, see

    Teleman J et al, "Numerical compression schemes for proteomics mass spectrometry data."
    Mol Cell Proteomics. 2014 Jun;13(6):1537-42. doi: 10.1074/mcp.O114.037879.

  */
  class OPENMS_DLLAPI MSNumpressCoder
  {

public:

    /// Names of compression schemes
    enum NumpressCompression {
      NONE, ///< No compression is applied
      LINEAR, ///< Linear (MS:1002312, MS-Numpress linear prediction compression)
      PIC, ///< Pic (MS:1002313, MS-Numpress positive integer compression)
      SLOF, ///< Slof (MS:1002314, MS-Numpress short logged float compression)
      SIZE_OF_NUMPRESSCOMPRESSION
    };
    static const std::string NamesOfNumpressCompression[SIZE_OF_NUMPRESSCOMPRESSION];

    /**
      @brief Configuration class for MSNumpress

      Contains configuration options for ms numpress
    */
    struct OPENMS_DLLAPI NumpressConfig
    {
      /**
        @brief fixed point for numpress algorithms

        Determines the accuracy of the encoding, is automatically estimated
        when estimate_fixed_point is set (only change this if you know what you
        are doing).

      */
      double numpressFixedPoint;
      /**
        @brief Check error tolerance after encoding

        Check error tolerance after encoding to ensure that the maximum error
        is abs(1.0-(encoded/decoded)) <= eps which is set here. In case it is
        set to 0, checking the encoding error is disabled. Note that this will
        slow down encoding substantially as all data needs to be encoded first
        and then decoded again.
      */
      double numpressErrorTolerance;
      /**
        @brief Which compression schema to use

        This is of type NumpressCompression  (see there)
      */
      NumpressCompression np_compression;
      /**
        @brief Whether to estimate the fixed point used for encoding (highly recommended)

        The fixed point determines the accuracy of the encoding and is
        automatically estimated when estimate_fixed_point is set to true.

        @note: only change this if you know what you are doing
      */
      bool estimate_fixed_point;
      /**
        @brief Desired mass accuracy for *linear* encoding

        This setting has no effect if set to -1, for example use 0.0001 for 0.2
        ppm accuracy @ 500 m/z. Does not affect other encoding schemes (pic or
        slof).
      */
      double linear_fp_mass_acc;

      NumpressConfig () :
        numpressFixedPoint(0.0),
        numpressErrorTolerance(BinaryDataEncoder_default_numpressErrorTolerance),
        np_compression(NONE),
        estimate_fixed_point(true),
        linear_fp_mass_acc(-1)
      {
      }

      /**
        @brief Set compression using a string mapping to enum NumpressCompression.

        @param compression A string from NamesOfNumpressCompression[]. Valid strings are "none", "linear", "pic" and "slof".

        @throws Exception::InvalidParameter if compression is unknown.
      */
      void setCompression(const std::string& compression)
      {
        const std::string* match = std::find(NamesOfNumpressCompression,
                                             NamesOfNumpressCompression + SIZE_OF_NUMPRESSCOMPRESSION, compression);

        if (match == NamesOfNumpressCompression + SIZE_OF_NUMPRESSCOMPRESSION) // == end()
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Value '" + compression + "' is not a valid Numpress compression scheme.");
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
    void encodeNP(const std::vector<double> & in,
                  String & result,
                  bool zlib_compression,
                  const NumpressConfig & config);

    /// encodeNP from a float (convert first to double)
    void encodeNP(const std::vector<float> & in,
                  String & result,
                  bool zlib_compression,
                  const NumpressConfig & config);

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
    void decodeNP(const String & in,
                  std::vector<double> & out,
                  bool zlib_compression,
                  const NumpressConfig & config);

    /**
     * @brief Encode the data vector "in" to a raw byte array
     *
     * @note In case of error, "result" is given back unmodified
     * @note The result is not a string but a raw byte array and may contain zero bytes
     *
     * This performs the raw numpress encoding on a set of data and does no
     * Base64 encoding on the result. Therefore the result string is likely
     * *unsafe* to handle and is a raw byte array.
     *
     * Please use the safe versions above unless you need access to the raw
     * byte arrays.
     *
     * @param in The vector of floating point numbers to be encoded
     * @param result The resulting string
     * @param config The numpress configuration defining the compression strategy
     *
    */
    void encodeNPRaw(const std::vector<double> & in,
                     String & result,
                     const NumpressConfig & config);

    /**
     * @brief Decode the raw byte array "in" to the result vector "out"
     *
     * @note The string in should *only* contain the data and _no_ extra
     * null terminating byte.
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
    void decodeNPRaw(const std::string & in,
                     std::vector<double> & out,
                     const NumpressConfig & config);

private:

    void decodeNPInternal_(const unsigned char* in, size_t in_size, std::vector<double>& out, const NumpressConfig & config);
  };

} //namespace OpenMS


