// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>

namespace OpenMS
{

  /**
    @brief Compresses and uncompresses byte buffers using Zstandard (zstd), optionally
           combined with the byte-shuffling and dictionary transforms defined for mzML.

    Static utility class implementing the zstd based binary data array compression
    methods recommended by the mzML specification (PSI-MS CV terms MS:1003780 to MS:1003785):

    - <b>zstd compression</b> (MS:1003780): the little-endian bytes of the array are compressed with zstd.
    - <b>byte-shuffled zstd compression</b> (MS:1003781): the bytes of the array elements are
      transposed ("shuffled") so that the i-th byte of every element is stored contiguously,
      followed by zstd compression. This usually compresses sorted data (e.g. m/z arrays) much better.
    - <b>dictionary-encoded zstd compression</b> (MS:1003782): the array is replaced by a sorted
      dictionary of its distinct values and an array of indices into that dictionary. Values and
      indices are byte-shuffled separately, then everything is compressed with zstd. This is well
      suited for arrays with many repeated values (e.g. ion mobility, charge states).
    - The MS-Numpress variants (MS:1003783 to MS:1003785) apply plain zstd compression to the
      output of the respective MS-Numpress encoder.

    The layout of a dictionary-encoded buffer is: an unsigned 64 bit little-endian integer holding
    the byte offset of the index array, an unsigned 64 bit little-endian integer holding the number
    of distinct values @em n, the byte-shuffled sorted distinct values, and finally the byte-shuffled
    indices. The indices use the smallest unsigned integer type able to address all values
    (8 bit if @em n < 2^8, 16 bit if @em n < 2^16, 32 bit if @em n < 2^32, otherwise 64 bit).

    All array transforms operate on raw byte buffers holding a contiguous array of fixed-size
    elements in little-endian byte order (the byte order mandated by mzML). The data is treated
    as raw bytes and may contain embedded zeros.

    @see https://github.com/mobiusklein/mzd.cpp for the reference implementation.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ZstdCompression
  {
public:

    /// Transform applied to the (little-endian) array bytes prior to zstd compression
    enum class ByteTransform
    {
      NONE,         ///< no transform, plain zstd compression (MS:1003780)
      BYTE_SHUFFLE, ///< byte shuffling followed by zstd compression (MS:1003781)
      DICTIONARY    ///< byte-shuffled dictionary encoding followed by zstd compression (MS:1003782)
    };

    /// Default zstd compression level (identical to zstd's own default)
    static constexpr int DEFAULT_LEVEL = 3;

    /**
      @brief Compress the @p in_length bytes pointed to by @p raw_data into @p compressed_data using zstd.

      A single zstd frame is written which records the uncompressed size in its header.

      @param[in]  raw_data        Pointer to the bytes to compress.
      @param[in]  in_length       Length of @p raw_data in bytes.
      @param[out] compressed_data Receives the compressed payload; any previous contents are replaced.
      @param[in]  level           The zstd compression level.

      @throws Exception::ConversionError if zstd reports a failure during compression.
    */
    static void compressData(const void* raw_data, size_t in_length, std::string& compressed_data, int level = DEFAULT_LEVEL);

    /**
      @brief Uncompress the zstd-compressed @p compressed_data.

      Multiple concatenated frames and frames without a recorded content size are supported.
      An empty input yields an empty output.

      @param[in]  compressed_data Pointer to the zstd-compressed bytes.
      @param[in]  nr_bytes        Length of @p compressed_data in bytes.
      @param[out] out             Receives the decompressed bytes; any previous contents are replaced.

      @throws Exception::ConversionError if the data is not valid zstd data or is truncated.
    */
    static void uncompressData(const void* compressed_data, size_t nr_bytes, std::string& out);

    /**
      @brief Byte-shuffle an array of @p element_size byte elements.

      The i-th byte of the j-th element is moved to position i * n + j (with @em n being the number of elements).

      @param[in]  data         Pointer to the array bytes.
      @param[in]  nr_bytes     Length of @p data in bytes; must be a multiple of @p element_size.
      @param[in]  element_size Size of a single array element in bytes.
      @param[out] out          Receives the shuffled bytes; any previous contents are replaced.

      @throws Exception::InvalidValue    if @p element_size is zero.
      @throws Exception::ConversionError if @p nr_bytes is not a multiple of @p element_size.
    */
    static void byteShuffle(const void* data, size_t nr_bytes, size_t element_size, std::string& out);

    /**
      @brief Reverse the byte shuffling done by byteShuffle().

      @param[in]  data         Pointer to the shuffled bytes.
      @param[in]  nr_bytes     Length of @p data in bytes; must be a multiple of @p element_size.
      @param[in]  element_size Size of a single array element in bytes.
      @param[out] out          Receives the restored array bytes; any previous contents are replaced.

      @throws Exception::InvalidValue    if @p element_size is zero.
      @throws Exception::ConversionError if @p nr_bytes is not a multiple of @p element_size.
    */
    static void byteUnshuffle(const void* data, size_t nr_bytes, size_t element_size, std::string& out);

    /**
      @brief Dictionary-encode an array of little-endian @p element_size byte elements (values and indices byte-shuffled).

      See the class documentation for the layout of the resulting buffer.

      @param[in]  data         Pointer to the array bytes (little-endian elements).
      @param[in]  nr_bytes     Length of @p data in bytes; must be a multiple of @p element_size.
      @param[in]  element_size Size of a single array element in bytes (1, 2, 4 or 8).
      @param[out] out          Receives the dictionary-encoded bytes; any previous contents are replaced.

      @throws Exception::InvalidValue    if @p element_size is not 1, 2, 4 or 8.
      @throws Exception::ConversionError if @p nr_bytes is not a multiple of @p element_size.
    */
    static void dictionaryEncode(const void* data, size_t nr_bytes, size_t element_size, std::string& out);

    /**
      @brief Decode a buffer created by dictionaryEncode() back into an array of little-endian @p element_size byte elements.

      An empty input yields an empty output.

      @param[in]  data         Pointer to the dictionary-encoded bytes.
      @param[in]  nr_bytes     Length of @p data in bytes.
      @param[in]  element_size Size of a single array element in bytes (1, 2, 4 or 8).
      @param[out] out          Receives the decoded array bytes; any previous contents are replaced.

      @throws Exception::InvalidValue    if @p element_size is not 1, 2, 4 or 8.
      @throws Exception::ConversionError if the buffer is malformed or its values do not have @p element_size bytes.
    */
    static void dictionaryDecode(const void* data, size_t nr_bytes, size_t element_size, std::string& out);

    /**
      @brief Apply @p transform to an array of little-endian @p element_size byte elements and compress the result with zstd.

      An empty input yields an empty output.

      @param[in]  data         Pointer to the array bytes (little-endian elements).
      @param[in]  nr_bytes     Length of @p data in bytes; must be a multiple of @p element_size.
      @param[in]  transform    The transform applied before compression.
      @param[in]  element_size Size of a single array element in bytes (ignored for ByteTransform::NONE).
      @param[out] out          Receives the compressed bytes; any previous contents are replaced.
      @param[in]  level        The zstd compression level.

      @throws Exception::InvalidValue    if @p element_size is invalid for @p transform.
      @throws Exception::ConversionError if the input size does not match @p element_size or compression fails.
    */
    static void encode(const void* data, size_t nr_bytes, ByteTransform transform, size_t element_size, std::string& out, int level = DEFAULT_LEVEL);

    /**
      @brief Uncompress zstd data and reverse @p transform, yielding an array of little-endian @p element_size byte elements.

      An empty input yields an empty output.

      @param[in]  data         Pointer to the compressed bytes.
      @param[in]  nr_bytes     Length of @p data in bytes.
      @param[in]  transform    The transform that was applied before compression.
      @param[in]  element_size Size of a single array element in bytes (ignored for ByteTransform::NONE).
      @param[out] out          Receives the decoded array bytes; any previous contents are replaced.

      @throws Exception::InvalidValue    if @p element_size is invalid for @p transform.
      @throws Exception::ConversionError if the data is malformed.
    */
    static void decode(const void* data, size_t nr_bytes, ByteTransform transform, size_t element_size, std::string& out);
  };

} // namespace OpenMS

