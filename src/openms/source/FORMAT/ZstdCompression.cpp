// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ZstdCompression.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <vector>

#include <zstd.h>

namespace OpenMS
{

  namespace
  {
    // Endianness-independent access to little-endian unsigned integers of arbitrary width
    template <typename UIntT>
    UIntT loadLittleEndian(const unsigned char* p)
    {
      UIntT value = 0;
      for (size_t b = 0; b < sizeof(UIntT); ++b)
      {
        value |= static_cast<UIntT>(p[b]) << (8 * b);
      }
      return value;
    }

    template <typename UIntT>
    void storeLittleEndian(UIntT value, unsigned char* p)
    {
      for (size_t b = 0; b < sizeof(UIntT); ++b)
      {
        p[b] = static_cast<unsigned char>(value >> (8 * b));
      }
    }

    void checkElementSize(size_t nr_bytes, size_t element_size, bool dictionary)
    {
      if (element_size == 0 ||
          (dictionary && element_size != 1 && element_size != 2 && element_size != 4 && element_size != 8))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported element size for zstd array transform", std::to_string(element_size));
      }
      if (nr_bytes % element_size != 0)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Array of " + std::to_string(nr_bytes) + " bytes is not a multiple of the element size " + std::to_string(element_size) + ".");
      }
    }

    /// Byte width of the dictionary indices for @p n_values distinct values (as defined in the mzML specification)
    size_t indexWidth(uint64_t n_values)
    {
      if (n_values < (uint64_t(1) << 8)) return 1;
      if (n_values < (uint64_t(1) << 16)) return 2;
      if (n_values < (uint64_t(1) << 32)) return 4;
      return 8;
    }

    uint64_t loadIndex(const unsigned char* p, size_t width)
    {
      switch (width)
      {
        case 1: return p[0];
        case 2: return loadLittleEndian<uint16_t>(p);
        case 4: return loadLittleEndian<uint32_t>(p);
        default: return loadLittleEndian<uint64_t>(p);
      }
    }

    void storeIndex(uint64_t index, unsigned char* p, size_t width)
    {
      switch (width)
      {
        case 1: p[0] = static_cast<unsigned char>(index); break;
        case 2: storeLittleEndian<uint16_t>(static_cast<uint16_t>(index), p); break;
        case 4: storeLittleEndian<uint32_t>(static_cast<uint32_t>(index), p); break;
        default: storeLittleEndian<uint64_t>(index, p); break;
      }
    }

    template <typename UIntT>
    void dictionaryEncodeImpl(const unsigned char* in, size_t count, std::string& out)
    {
      // interpret the elements as unsigned integers (by their bit pattern) and collect the sorted distinct values
      std::vector<UIntT> values(count);
      for (size_t j = 0; j < count; ++j)
      {
        values[j] = loadLittleEndian<UIntT>(in + j * sizeof(UIntT));
      }
      std::vector<UIntT> dictionary(values);
      std::sort(dictionary.begin(), dictionary.end());
      dictionary.erase(std::unique(dictionary.begin(), dictionary.end()), dictionary.end());

      const uint64_t n_values = dictionary.size();
      const size_t width = indexWidth(n_values);
      const size_t values_bytes = dictionary.size() * sizeof(UIntT);
      const uint64_t offset = 16 + values_bytes;

      // little-endian dictionary values and indices (shuffled below)
      std::string dict_bytes(values_bytes, '\0');
      for (size_t i = 0; i < dictionary.size(); ++i)
      {
        storeLittleEndian<UIntT>(dictionary[i], reinterpret_cast<unsigned char*>(&dict_bytes[i * sizeof(UIntT)]));
      }
      std::string index_bytes(count * width, '\0');
      for (size_t j = 0; j < count; ++j)
      {
        const uint64_t index = std::lower_bound(dictionary.begin(), dictionary.end(), values[j]) - dictionary.begin();
        storeIndex(index, reinterpret_cast<unsigned char*>(&index_bytes[j * width]), width);
      }

      std::string shuffled_values, shuffled_indices;
      ZstdCompression::byteShuffle(dict_bytes.data(), dict_bytes.size(), sizeof(UIntT), shuffled_values);
      ZstdCompression::byteShuffle(index_bytes.data(), index_bytes.size(), width, shuffled_indices);

      out.resize(16);
      storeLittleEndian<uint64_t>(offset, reinterpret_cast<unsigned char*>(&out[0]));
      storeLittleEndian<uint64_t>(n_values, reinterpret_cast<unsigned char*>(&out[8]));
      out.append(shuffled_values);
      out.append(shuffled_indices);
    }

    [[noreturn]] void throwZstdError(size_t code, const std::string& what)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       what + ": " + ZSTD_getErrorName(code));
    }

    struct DCtxDeleter
    {
      void operator()(ZSTD_DCtx* ctx) const { ZSTD_freeDCtx(ctx); }
    };
  }

  void ZstdCompression::compressData(const void* raw_data, size_t in_length, std::string& compressed_data, int level)
  {
    compressed_data.clear();
    const size_t bound = ZSTD_compressBound(in_length);
    if (ZSTD_isError(bound))
    {
      throwZstdError(bound, "zstd compression failed");
    }
    compressed_data.resize(bound);
    const size_t used = ZSTD_compress(&compressed_data[0], bound, raw_data, in_length, level);
    if (ZSTD_isError(used))
    {
      throwZstdError(used, "zstd compression failed");
    }
    compressed_data.resize(used);
  }

  void ZstdCompression::uncompressData(const void* compressed_data, size_t nr_bytes, std::string& out)
  {
    out.clear();
    if (nr_bytes == 0)
    {
      return;
    }

    // Use the content size stored in the (first) frame header as initial buffer size. Streaming
    // decompression is used so that frames without a recorded size, multiple frames and headers
    // claiming an incorrect size are handled gracefully: the buffer simply grows when needed.
    size_t capacity = ZSTD_DStreamOutSize();
    const unsigned long long content_size = ZSTD_getFrameContentSize(compressed_data, nr_bytes);
    if (content_size == ZSTD_CONTENTSIZE_ERROR)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "zstd decompression failed: input is not valid zstd data");
    }
    if (content_size != ZSTD_CONTENTSIZE_UNKNOWN)
    {
      // guard against absurd sizes in (malformed) frame headers; the buffer grows if required
      capacity = static_cast<size_t>(std::min<unsigned long long>(std::max<unsigned long long>(content_size, 1), 1ull << 28));
    }

    std::unique_ptr<ZSTD_DCtx, DCtxDeleter> dctx(ZSTD_createDCtx());
    if (!dctx)
    {
      throw Exception::OutOfMemory(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sizeof(ZSTD_DCtx*));
    }

    out.resize(capacity);
    size_t produced = 0;
    ZSTD_inBuffer input = {compressed_data, nr_bytes, 0};
    while (true)
    {
      if (produced == out.size())
      {
        out.resize(out.size() * 2);
      }
      ZSTD_outBuffer output = {&out[produced], out.size() - produced, 0};
      const size_t ret = ZSTD_decompressStream(dctx.get(), &output, &input);
      if (ZSTD_isError(ret))
      {
        throwZstdError(ret, "zstd decompression failed");
      }
      produced += output.pos;
      if (input.pos == input.size)
      {
        if (ret == 0)
        {
          break; // all frames completely decoded
        }
        if (output.pos < output.size)
        {
          // no more input, but the decoder expects more data and has nothing left to flush
          throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "zstd decompression failed: truncated input");
        }
      }
    }
    out.resize(produced);
  }

  void ZstdCompression::byteShuffle(const void* data, size_t nr_bytes, size_t element_size, std::string& out)
  {
    checkElementSize(nr_bytes, element_size, false);
    const unsigned char* in = static_cast<const unsigned char*>(data);
    const size_t count = nr_bytes / element_size;
    out.resize(nr_bytes);
    for (size_t i = 0; i < element_size; ++i)
    {
      char* dest = &out[0] + i * count;
      for (size_t j = 0; j < count; ++j)
      {
        dest[j] = static_cast<char>(in[j * element_size + i]);
      }
    }
  }

  void ZstdCompression::byteUnshuffle(const void* data, size_t nr_bytes, size_t element_size, std::string& out)
  {
    checkElementSize(nr_bytes, element_size, false);
    const unsigned char* in = static_cast<const unsigned char*>(data);
    const size_t count = nr_bytes / element_size;
    out.resize(nr_bytes);
    for (size_t i = 0; i < element_size; ++i)
    {
      const unsigned char* src = in + i * count;
      for (size_t j = 0; j < count; ++j)
      {
        out[j * element_size + i] = static_cast<char>(src[j]);
      }
    }
  }

  void ZstdCompression::dictionaryEncode(const void* data, size_t nr_bytes, size_t element_size, std::string& out)
  {
    checkElementSize(nr_bytes, element_size, true);
    out.clear();
    const unsigned char* in = static_cast<const unsigned char*>(data);
    const size_t count = nr_bytes / element_size;
    switch (element_size)
    {
      case 1: dictionaryEncodeImpl<uint8_t>(in, count, out); break;
      case 2: dictionaryEncodeImpl<uint16_t>(in, count, out); break;
      case 4: dictionaryEncodeImpl<uint32_t>(in, count, out); break;
      default: dictionaryEncodeImpl<uint64_t>(in, count, out); break;
    }
  }

  void ZstdCompression::dictionaryDecode(const void* data, size_t nr_bytes, size_t element_size, std::string& out)
  {
    checkElementSize(0, element_size, true);
    out.clear();
    if (nr_bytes == 0)
    {
      return;
    }
    if (nr_bytes < 16)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Malformed dictionary-encoded array: buffer of " + std::to_string(nr_bytes) + " bytes is shorter than its 16 byte header.");
    }
    const unsigned char* in = static_cast<const unsigned char*>(data);
    const uint64_t offset = loadLittleEndian<uint64_t>(in);
    const uint64_t n_values = loadLittleEndian<uint64_t>(in + 8);
    if (offset < 16 || offset > nr_bytes)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Malformed dictionary-encoded array: index offset " + std::to_string(offset) + " is outside of the buffer of " + std::to_string(nr_bytes) + " bytes.");
    }
    const uint64_t values_bytes = offset - 16;
    if (n_values == 0)
    {
      if (values_bytes != 0 || offset != nr_bytes)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Malformed dictionary-encoded array: empty dictionary, but buffer contains data.");
      }
      return;
    }
    if (values_bytes % element_size != 0 || values_bytes / element_size != n_values)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Malformed dictionary-encoded array: " + std::to_string(n_values) + " dictionary values do not occupy " +
                                       std::to_string(values_bytes) + " bytes with an element size of " + std::to_string(element_size) + " bytes.");
    }
    const size_t width = indexWidth(n_values);
    const size_t index_bytes = nr_bytes - static_cast<size_t>(offset);
    if (index_bytes % width != 0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Malformed dictionary-encoded array: index region of " + std::to_string(index_bytes) + " bytes is not a multiple of the index width " + std::to_string(width) + ".");
    }

    std::string dictionary, indices;
    byteUnshuffle(in + 16, static_cast<size_t>(values_bytes), element_size, dictionary);
    byteUnshuffle(in + offset, index_bytes, width, indices);

    const size_t count = index_bytes / width;
    out.resize(count * element_size);
    const unsigned char* idx_ptr = reinterpret_cast<const unsigned char*>(indices.data());
    for (size_t j = 0; j < count; ++j)
    {
      const uint64_t index = loadIndex(idx_ptr + j * width, width);
      if (index >= n_values)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Malformed dictionary-encoded array: index " + std::to_string(index) + " exceeds the dictionary size " + std::to_string(n_values) + ".");
      }
      std::memcpy(&out[j * element_size], &dictionary[static_cast<size_t>(index) * element_size], element_size);
    }
  }

  void ZstdCompression::encode(const void* data, size_t nr_bytes, ByteTransform transform, size_t element_size, std::string& out, int level)
  {
    out.clear();
    if (nr_bytes == 0)
    {
      return;
    }
    std::string transformed;
    switch (transform)
    {
      case ByteTransform::NONE:
        compressData(data, nr_bytes, out, level);
        return;
      case ByteTransform::BYTE_SHUFFLE:
        byteShuffle(data, nr_bytes, element_size, transformed);
        break;
      case ByteTransform::DICTIONARY:
        dictionaryEncode(data, nr_bytes, element_size, transformed);
        break;
    }
    compressData(transformed.data(), transformed.size(), out, level);
  }

  void ZstdCompression::decode(const void* data, size_t nr_bytes, ByteTransform transform, size_t element_size, std::string& out)
  {
    out.clear();
    if (nr_bytes == 0)
    {
      return;
    }
    if (transform == ByteTransform::NONE)
    {
      uncompressData(data, nr_bytes, out);
      return;
    }
    std::string uncompressed;
    uncompressData(data, nr_bytes, uncompressed);
    if (transform == ByteTransform::BYTE_SHUFFLE)
    {
      byteUnshuffle(uncompressed.data(), uncompressed.size(), element_size, out);
    }
    else
    {
      dictionaryDecode(uncompressed.data(), uncompressed.size(), element_size, out);
    }
  }

} // namespace OpenMS
