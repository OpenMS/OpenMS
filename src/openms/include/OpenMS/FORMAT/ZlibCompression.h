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

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <vector>
#include <QByteArray>

namespace OpenMS
{
  class String;

  /**
    * @brief Compresses and uncompresses data using zlib
    *
    * @note The 'strings' here are not really null-terminated but rather
    * containers of data. If you want safe conversions, use QtByteArray.
    * 
  */
  class OPENMS_DLLAPI ZlibCompression
  {
public:

    /**
      * @brief Compresses data using zlib directly
      *
      * @param raw_data Data to be compressed
      * @param compressed_data Compressed result data
      * 
    */
    static void compressString(std::string& raw_data, std::string& compressed_data);

    /**
      * @brief Compresses data using Qt
      *
      * @param raw_data Data to be compressed
      * @param compressed_data Compressed result data
      * 
    */
    static void compressString(const QByteArray& raw_data, QByteArray& compressed_data);

    /**
      * @brief Uncompresses data using Qt (wrapper around Qt function)
      *
      * @param compressed_data Compressed data
      * @param nr_bytes Number of bytes in compressed data
      * @param raw_data Uncompressed result data
      * 
    */
    static void uncompressString(const void * compressed_data, size_t nr_bytes, std::string& raw_data);

    /**
      * @brief Uncompresses data using Qt
      *
      * @param compressed_data Compressed data
      * @param raw_data Uncompressed result data
      * 
    */
    static void uncompressString(const QByteArray& compressed_data, QByteArray& raw_data);

  };

} // namespace OpenMS


