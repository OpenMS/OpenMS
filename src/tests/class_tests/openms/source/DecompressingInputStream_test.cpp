// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/DecompressingInputStream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/FORMAT/GzipInputStream.h>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>

using namespace OpenMS;

namespace
{
  // An in-memory stream that hands out at most chunk bytes per read, to exercise buffer boundaries
  class MemoryStream : public InputStream
  {
public:
    MemoryStream(std::string data, Size chunk) : data_(std::move(data)), chunk_(chunk) {}
    UInt64 curPos() const override { return position_; }
    Size readBytes(unsigned char* to_fill, Size max_to_read) override
    {
      const Size count = std::min({max_to_read, chunk_, data_.size() - position_});
      std::copy_n(data_.data() + position_, count, to_fill);
      position_ += count;
      return count;
    }
private:
    std::string data_;
    Size chunk_;
    Size position_ = 0;
  };

  std::string readFile(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }

  std::string readAll(InputStream& stream, Size block = 4096)
  {
    std::string result;
    std::string buffer(block, '\0');
    while (const Size count = stream.readBytes(reinterpret_cast<unsigned char*>(buffer.data()), block))
    {
      result.append(buffer.data(), count);
    }
    return result;
  }

  std::string decompress(const std::string& data, DecompressingInputStream::Format format, Size chunk = 7, Size block = 4096)
  {
    DecompressingInputStream stream(std::make_unique<MemoryStream>(data, chunk), format, "test");
    return readAll(stream, block);
  }
}

START_TEST(DecompressingInputStream, "$Id$")

const std::string gz_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz");
const std::string bz2_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2");
const std::string gz = readFile(gz_path);
const std::string bz2 = readFile(bz2_path);
GzipInputStream gzip_reference(gz_path);
const std::string gz_plain = readAll(gzip_reference);
Bzip2InputStream bzip2_reference(bz2_path);
const std::string bz2_plain = readAll(bzip2_reference);
using Format = DecompressingInputStream::Format;

START_SECTION((static std::optional<Format> detect(const unsigned char* head, Size size)))
{
  auto detect = [](const std::string& head)
  {
    return DecompressingInputStream::detect(reinterpret_cast<const unsigned char*>(head.data()), head.size());
  };
  TEST_TRUE(detect(gz) == Format::GZIP)
  TEST_TRUE(detect(bz2) == Format::BZIP2)
  TEST_TRUE(detect("<?xml") == std::nullopt)
  TEST_TRUE(detect("BZ") == std::nullopt) // too short for bzip2
  TEST_TRUE(detect(std::string("\x1f", 1)) == std::nullopt)
  TEST_TRUE(detect("PK\x03\x04") == std::nullopt) // ZIP is not a stream format
}
END_SECTION

START_SECTION((DecompressingInputStream(std::unique_ptr<InputStream> source, Format format, const std::string& name)))
{
  TEST_EQUAL(gz_plain.empty(), false)
  TEST_EQUAL(decompress(gz, Format::GZIP), gz_plain)
  TEST_EQUAL(decompress(bz2, Format::BZIP2), bz2_plain)
}
END_SECTION

START_SECTION((Size readBytes(unsigned char* to_fill, Size max_to_read)))
{
  // whole-buffer and byte-by-byte reads
  TEST_EQUAL(decompress(gz, Format::GZIP, 1 << 20, 1 << 20), gz_plain)
  TEST_EQUAL(decompress(gz, Format::GZIP, 1, 1), gz_plain)
  TEST_EQUAL(decompress(bz2, Format::BZIP2, 3, 5), bz2_plain)

  // concatenated members (pigz, pbzip2) are read in sequence
  TEST_EQUAL(decompress(gz + gz, Format::GZIP), gz_plain + gz_plain)
  TEST_EQUAL(decompress(bz2 + bz2 + bz2, Format::BZIP2), bz2_plain + bz2_plain + bz2_plain)

  // other trailing bytes are ignored, as gzread() does
  TEST_EQUAL(decompress(gz + std::string(5, '\0'), Format::GZIP), gz_plain)
  TEST_EQUAL(decompress(bz2 + "x", Format::BZIP2), bz2_plain)

  // truncated or corrupt data is an error, not a silently shorter document
  TEST_EXCEPTION(Exception::ParseError, decompress(gz.substr(0, gz.size() / 2), Format::GZIP))
  TEST_EXCEPTION(Exception::ParseError, decompress(bz2.substr(0, bz2.size() - 10), Format::BZIP2))
  std::string corrupt = gz;
  for (Size i = corrupt.size() / 2; i < corrupt.size() / 2 + 40; ++i) { corrupt[i] = static_cast<char>(~corrupt[i]); }
  TEST_EXCEPTION(Exception::ParseError, decompress(corrupt, Format::GZIP))
  TEST_EXCEPTION(Exception::ParseError, decompress("<?xml version", Format::BZIP2))

  // nothing to read
  DecompressingInputStream stream(std::make_unique<MemoryStream>(gz, 7), Format::GZIP, "test");
  unsigned char byte = 0;
  TEST_EQUAL(stream.readBytes(&byte, 0), 0)
}
END_SECTION

START_SECTION((UInt64 curPos() const))
{
  DecompressingInputStream stream(std::make_unique<MemoryStream>(gz, 7), Format::GZIP, "test");
  TEST_EQUAL(stream.curPos(), 0)
  unsigned char buffer[100];
  TEST_EQUAL(stream.readBytes(buffer, sizeof(buffer)), sizeof(buffer))
  TEST_EQUAL(stream.curPos(), sizeof(buffer))
  readAll(stream);
  TEST_EQUAL(stream.curPos(), gz_plain.size())
}
END_SECTION

END_TEST
