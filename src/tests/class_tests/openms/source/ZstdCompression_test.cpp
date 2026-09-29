// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ZstdCompression.h>
///////////////////////////

#include <cstdint>
#include <cstring>
#include <vector>

using namespace OpenMS;

namespace
{
  using BT = ZstdCompression::ByteTransform;

  // test vectors generated independently with Python (numpy + zstandard) following the
  // reference code of the PSI mzML binary array compression specification
  // double array {100.5, 200.25, 300.125, 100.5, 400.0} and int32 array {1, 2, 3, -4, 1, 2}
  const std::string d_plain(
    "\x28\xb5\x2f\xfd\x20\x28\xe5\x00\x00\x80\x00\x00\x00\x00\x00\x20\x59\x40\x00\x08\x69\xc2\x72\x00"
    "\x79\x40\x04\x00\x3b\x09\x00\x03\xaf\x0c\x0e\xc0\x11", 37);
  const std::string d_shuffle(
    "\x28\xb5\x2f\xfd\x20\x28\xbd\x00\x00\x88\x00\x00\x20\x08\xc2\x20\x00\x59\x69\x72\x59\x79\x40\x40"
    "\x40\x40\x40\x01\x00\x35\xc0\x02", 32);
  const std::string d_dict(
    "\x28\xb5\x2f\xfd\x20\x35\xed\x00\x00\xa8\x30\x00\x04\x00\x20\x08\xc2\x00\x59\x69\x72\x79\x40\x40"
    "\x40\x40\x00\x01\x02\x00\x03\x02\x00\xc0\x62\x0c\x60\x01", 38);
  const std::string i_dict(
    "\x28\xb5\x2f\xfd\x20\x26\xdd\x00\x00\x88\x20\x00\x04\x00\x01\x02\x03\xfc\x00\x00\x00\xff\x01\x02"
    "\x03\x00\x01\x03\x00\xbb\xc2\xc1\xc0\x03\xc0\x02", 36);
  const std::string i_dict_raw(
    "\x20\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x01\x02\x03\xfc\x00\x00\x00\xff"
    "\x00\x00\x00\xff\x00\x00\x00\xff\x00\x01\x02\x03\x00\x01", 38);
  const std::string d_shuffle_raw(
    "\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    "\x00\x20\x08\xc2\x20\x00\x59\x69\x72\x59\x79\x40\x40\x40\x40\x40", 40);
  const std::string d_multi(
    "\x28\xb5\x2f\xfd\x00\x00\x81\x00\x00\x00\x00\x00\x00\x00\x20\x59\x40\x00\x00\x00\x00\x00\x08\x69"
    "\x40\x28\xb5\x2f\xfd\x20\x18\xc1\x00\x00\x00\x00\x00\x00\x00\xc2\x72\x40\x00\x00\x00\x00\x00\x20"
    "\x59\x40\x00\x00\x00\x00\x00\x00\x79\x40", 58);

  template <typename T>
  std::vector<T> toVector(const std::string& bytes)
  {
    std::vector<T> result(bytes.size() / sizeof(T));
    if (!result.empty())
    {
      std::memcpy(result.data(), bytes.data(), result.size() * sizeof(T));
    }
    return result;
  }

  template <typename T>
  std::string toBytes(const std::vector<T>& values)
  {
    return std::string(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(T));
  }

  // Hand-built dictionary-encoded buffer (MS:1003782 layout) of 16 bit values with @p index_width byte
  // indices: n_values distinct values (0, 3, 6, ...) and count elements referring to value (j % n_values).
  // Returns the buffer and stores the plain little-endian array in @p plain.
  std::string buildDictionary(size_t n_values, size_t count, size_t index_width, std::string& plain)
  {
    std::string values(n_values * 2, '\0');
    for (size_t i = 0; i < n_values; ++i)
    {
      values[2 * i] = static_cast<char>((i * 3) & 0xFF);
      values[2 * i + 1] = static_cast<char>((i * 3) >> 8);
    }
    std::string indices(count * index_width, '\0');
    plain.clear();
    for (size_t j = 0; j < count; ++j)
    {
      const size_t index = j % n_values;
      for (size_t b = 0; b < index_width; ++b)
      {
        indices[j * index_width + b] = static_cast<char>((index >> (8 * b)) & 0xFF);
      }
      plain.append(values, 2 * index, 2);
    }
    std::string shuffled_values, shuffled_indices;
    ZstdCompression::byteShuffle(values.data(), values.size(), 2, shuffled_values);
    ZstdCompression::byteShuffle(indices.data(), indices.size(), index_width, shuffled_indices);
    std::string header(16, '\0');
    const uint64_t offset = 16 + values.size();
    for (size_t b = 0; b < 8; ++b)
    {
      header[b] = static_cast<char>((offset >> (8 * b)) & 0xFF);
      header[8 + b] = static_cast<char>((uint64_t(n_values) >> (8 * b)) & 0xFF);
    }
    return header + shuffled_values + shuffled_indices;
  }
}

START_TEST(ZstdCompression, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// note: the test vectors assume a little-endian platform
const std::vector<double> d_values = {100.5, 200.25, 300.125, 100.5, 400.0};
const std::vector<Int32> i_values = {1, 2, 3, -4, 1, 2};

START_SECTION((static void compressData(const void* raw_data, size_t in_length, std::string& compressed_data, int level = DEFAULT_LEVEL)))
{
  const std::string raw = "Freude, schoner Gotterfunken, Tochter aus Elysium, Freude, schoner Gotterfunken, Tochter aus Elysium!";
  std::string compressed;
  ZstdCompression::compressData(raw.data(), raw.size(), compressed);
  TEST_TRUE(compressed.size() > 4)
  TEST_TRUE(compressed.size() < raw.size())
  // zstd magic number
  TEST_EQUAL(compressed.substr(0, 4), std::string("\x28\xb5\x2f\xfd", 4))
  std::string uncompressed;
  ZstdCompression::uncompressData(compressed.data(), compressed.size(), uncompressed);
  TEST_EQUAL(uncompressed, raw)

  // other compression levels
  ZstdCompression::compressData(raw.data(), raw.size(), compressed, 19);
  ZstdCompression::uncompressData(compressed.data(), compressed.size(), uncompressed);
  TEST_EQUAL(uncompressed, raw)

  // empty input
  ZstdCompression::compressData(raw.data(), 0, compressed);
  ZstdCompression::uncompressData(compressed.data(), compressed.size(), uncompressed);
  TEST_EQUAL(uncompressed.size(), 0)
}
END_SECTION

START_SECTION((static void uncompressData(const void* compressed_data, size_t nr_bytes, std::string& out, size_t expected_size = 0)))
{
  std::string out;
  ZstdCompression::uncompressData(d_plain.data(), d_plain.size(), out);
  TEST_EQUAL(out, toBytes(d_values))

  // multiple frames, the first one without recorded content size
  ZstdCompression::uncompressData(d_multi.data(), d_multi.size(), out);
  TEST_EQUAL(out, toBytes(d_values))

  // empty input
  std::string empty;
  ZstdCompression::uncompressData(empty.data(), 0, empty);
  TEST_EQUAL(empty.size(), 0)

  // large data (exceeding the streaming buffer size)
  std::string large;
  for (Size i = 0; i < 500000; ++i)
  {
    large += std::to_string(i % 977);
  }
  std::string compressed, uncompressed;
  ZstdCompression::compressData(large.data(), large.size(), compressed);
  ZstdCompression::uncompressData(compressed.data(), compressed.size(), uncompressed);
  TEST_EQUAL(uncompressed == large, true)

  // the expected size only limits the initial allocation: too small, exact or too large sizes yield the same result
  for (size_t expected_size : {size_t(1), size_t(1000), large.size(), 10 * large.size()})
  {
    ZstdCompression::uncompressData(compressed.data(), compressed.size(), uncompressed, expected_size);
    TEST_EQUAL(uncompressed == large, true)
  }

  // invalid and truncated data
  const std::string invalid = "this is not zstd";
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::uncompressData(invalid.data(), invalid.size(), uncompressed))
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::uncompressData(compressed.data(), compressed.size() / 2, uncompressed))
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::uncompressData(compressed.data(), compressed.size() / 2, uncompressed, large.size()))
}
END_SECTION

START_SECTION((static void byteShuffle(const void* data, size_t nr_bytes, size_t element_size, std::string& out)))
{
  const std::vector<UInt32> values = {0x04030201, 0x08070605, 0x0C0B0A09};
  std::string shuffled;
  ZstdCompression::byteShuffle(values.data(), values.size() * sizeof(UInt32), sizeof(UInt32), shuffled);
  TEST_EQUAL(shuffled, std::string("\x01\x05\x09\x02\x06\x0A\x03\x07\x0B\x04\x08\x0C", 12))

  ZstdCompression::byteShuffle(d_values.data(), d_values.size() * sizeof(double), sizeof(double), shuffled);
  TEST_EQUAL(shuffled, d_shuffle_raw)

  // element size 1 is the identity
  const std::string text = "abcdef";
  ZstdCompression::byteShuffle(text.data(), text.size(), 1, shuffled);
  TEST_EQUAL(shuffled, text)

  // arbitrary element sizes are supported
  ZstdCompression::byteShuffle(text.data(), text.size(), 3, shuffled);
  TEST_EQUAL(shuffled, "adbecf")

  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::byteShuffle(text.data(), text.size(), 4, shuffled))
  TEST_EXCEPTION(Exception::InvalidValue, ZstdCompression::byteShuffle(text.data(), text.size(), 0, shuffled))
}
END_SECTION

START_SECTION((static void byteUnshuffle(const void* data, size_t nr_bytes, size_t element_size, std::string& out)))
{
  const std::string shuffled("\x01\x05\x09\x02\x06\x0A\x03\x07\x0B\x04\x08\x0C", 12);
  std::string out;
  ZstdCompression::byteUnshuffle(shuffled.data(), shuffled.size(), 4, out);
  TEST_EQUAL(out, std::string("\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0A\x0B\x0C", 12))
  ZstdCompression::byteUnshuffle(shuffled.data(), shuffled.size(), 2, out);
  std::string back;
  ZstdCompression::byteShuffle(out.data(), out.size(), 2, back);
  TEST_EQUAL(back, shuffled)
  // 12 bytes are not a multiple of 8
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::byteUnshuffle(shuffled.data(), shuffled.size(), 8, out))
}
END_SECTION

START_SECTION((static void dictionaryEncode(const void* data, size_t nr_bytes, size_t element_size, std::string& out)))
{
  std::string encoded;
  ZstdCompression::dictionaryEncode(i_values.data(), i_values.size() * sizeof(Int32), sizeof(Int32), encoded);
  TEST_EQUAL(encoded, i_dict_raw)

  // the width of the indices depends on the number of unique values
  for (Size n : {1, 255, 256, 65535, 65536})
  {
    std::vector<double> values;
    for (Size k = 0; k < 2 * n; ++k)
    {
      values.push_back(0.25 * static_cast<double>(k % n));
    }
    ZstdCompression::dictionaryEncode(values.data(), values.size() * sizeof(double), sizeof(double), encoded);
    const Size index_width = n < 256 ? 1 : (n < 65536 ? 2 : 4);
    TEST_EQUAL(encoded.size(), 16 + n * sizeof(double) + 2 * n * index_width)
    std::string decoded;
    ZstdCompression::dictionaryDecode(encoded.data(), encoded.size(), sizeof(double), decoded);
    TEST_EQUAL(decoded == toBytes(values), true)
  }

  // empty input: header only
  const std::vector<double> no_values;
  ZstdCompression::dictionaryEncode(no_values.data(), 0, sizeof(double), encoded);
  TEST_EQUAL(encoded.size(), 16)
  std::string decoded;
  ZstdCompression::dictionaryDecode(encoded.data(), encoded.size(), sizeof(double), decoded);
  TEST_EQUAL(decoded.size(), 0)

  // only element sizes of 1, 2, 4 and 8 bytes are supported
  TEST_EXCEPTION(Exception::InvalidValue, ZstdCompression::dictionaryEncode(i_values.data(), 3, 3, encoded))
}
END_SECTION

START_SECTION((static void dictionaryDecode(const void* data, size_t nr_bytes, size_t element_size, std::string& out, size_t array_length = 0)))
{
  std::string decoded;
  ZstdCompression::dictionaryDecode(i_dict_raw.data(), i_dict_raw.size(), sizeof(Int32), decoded);
  TEST_EQUAL(decoded, toBytes(i_values))
  // malformed buffers
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::dictionaryDecode("\x00\x01", 2, sizeof(Int32), decoded))
  std::string wrong_offset = i_dict_raw;
  wrong_offset[0] = '\x7F';
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::dictionaryDecode(wrong_offset.data(), wrong_offset.size(), sizeof(Int32), decoded))
  // element size 8 does not match the dictionary of 4 values starting at offset 32
  TEST_EXCEPTION(Exception::ConversionError, ZstdCompression::dictionaryDecode(i_dict_raw.data(), i_dict_raw.size(), sizeof(double), decoded))

  // a (correct) array length does not change the result for buffers following the specification
  ZstdCompression::dictionaryDecode(i_dict_raw.data(), i_dict_raw.size(), sizeof(Int32), decoded, i_values.size());
  TEST_EQUAL(decoded, toBytes(i_values))

  // index widths at the type boundaries: the specification uses 8 bit indices for up to 255 values,
  // the Rust library mzdata (mobiusklein/mzdata) 8 bit indices for up to 256 values and mzd.cpp's
  // reader 16 bit indices for 255 values
  std::string plain;
  const size_t count = 300;
  // 256 values, 8 bit indices (Rust library mzdata, mobiusklein/mzdata): decodable with the array length
  std::string dict = buildDictionary(256, count, 1, plain);
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded, count);
  TEST_EQUAL(decoded == plain, true)
  // 256 values, 16 bit indices (specification): decodable with and without the array length
  dict = buildDictionary(256, count, 2, plain);
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded, count);
  TEST_EQUAL(decoded == plain, true)
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded);
  TEST_EQUAL(decoded == plain, true)
  // 255 values, 16 bit indices (mzd.cpp reader): decodable with the array length
  dict = buildDictionary(255, count, 2, plain);
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded, count);
  TEST_EQUAL(decoded == plain, true)
  // 255 values, 8 bit indices (specification): decodable with and without the array length
  dict = buildDictionary(255, count, 1, plain);
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded, count);
  TEST_EQUAL(decoded == plain, true)
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded);
  TEST_EQUAL(decoded == plain, true)
  // an array length that does not match the index region falls back to the specification's width
  ZstdCompression::dictionaryDecode(dict.data(), dict.size(), 2, decoded, count + 1);
  TEST_EQUAL(decoded == plain, true)
}
END_SECTION

START_SECTION((static void encode(const void* data, size_t nr_bytes, ByteTransform transform, size_t element_size, std::string& out, int level = DEFAULT_LEVEL)))
{
  std::vector<float> values;
  for (Size k = 0; k < 10000; ++k)
  {
    values.push_back(static_cast<float>(k % 100) * 1.5f);
  }
  for (BT transform : {BT::NONE, BT::BYTE_SHUFFLE, BT::DICTIONARY})
  {
    std::string encoded, decoded;
    ZstdCompression::encode(values.data(), values.size() * sizeof(float), transform, sizeof(float), encoded);
    TEST_TRUE(encoded.size() < values.size() * sizeof(float))
    ZstdCompression::decode(encoded.data(), encoded.size(), transform, sizeof(float), decoded);
    TEST_EQUAL(toVector<float>(decoded) == values, true)
  }
}
END_SECTION

START_SECTION((static void decode(const void* data, size_t nr_bytes, ByteTransform transform, size_t element_size, std::string& out, size_t array_length = 0)))
{
  std::string decoded;
  ZstdCompression::decode(d_plain.data(), d_plain.size(), BT::NONE, sizeof(double), decoded);
  TEST_EQUAL(toVector<double>(decoded) == d_values, true)
  ZstdCompression::decode(d_shuffle.data(), d_shuffle.size(), BT::BYTE_SHUFFLE, sizeof(double), decoded);
  TEST_EQUAL(toVector<double>(decoded) == d_values, true)
  ZstdCompression::decode(d_dict.data(), d_dict.size(), BT::DICTIONARY, sizeof(double), decoded);
  TEST_EQUAL(toVector<double>(decoded) == d_values, true)
  ZstdCompression::decode(i_dict.data(), i_dict.size(), BT::DICTIONARY, sizeof(Int32), decoded);
  TEST_EQUAL(toVector<Int32>(decoded) == i_values, true)

  // with array length
  ZstdCompression::decode(d_plain.data(), d_plain.size(), BT::NONE, sizeof(double), decoded, d_values.size());
  TEST_EQUAL(toVector<double>(decoded) == d_values, true)
  ZstdCompression::decode(d_shuffle.data(), d_shuffle.size(), BT::BYTE_SHUFFLE, sizeof(double), decoded, d_values.size());
  TEST_EQUAL(toVector<double>(decoded) == d_values, true)
  ZstdCompression::decode(i_dict.data(), i_dict.size(), BT::DICTIONARY, sizeof(Int32), decoded, i_values.size());
  TEST_EQUAL(toVector<Int32>(decoded) == i_values, true)

  // the array length is passed on to the dictionary decoder (256 values with 8 bit indices)
  std::string plain;
  const std::string dict = buildDictionary(256, 300, 1, plain);
  std::string compressed;
  ZstdCompression::compressData(dict.data(), dict.size(), compressed);
  ZstdCompression::decode(compressed.data(), compressed.size(), BT::DICTIONARY, 2, decoded, 300);
  TEST_EQUAL(decoded == plain, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
