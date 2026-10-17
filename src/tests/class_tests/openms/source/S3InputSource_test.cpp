// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Init.h>
#include "FORMAT/DATAACCESS/S3InputSource.h"

#include <arrow/io/memory.h>
#include <arrow/buffer.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <xercesc/util/XMLString.hpp>
#include <algorithm>
#include <array>

using namespace OpenMS;

namespace
{
  // Force partial reads even when Xerces (or the decompressor) asks for more.
  class ChunkedStream : public arrow::io::InputStream
  {
  public:
    explicit ChunkedStream(const std::string& data, int64_t fail_after = -1)
      : reader_(arrow::Buffer::FromString(data)), fail_after_(fail_after)
    {
    }

    arrow::Result<int64_t> Read(int64_t nbytes, void* out) override
    {
      auto position = reader_.Tell();
      if (!position.ok()) return position.status();
      if (fail_after_ >= 0 && *position >= fail_after_)
      {
        return arrow::Status::IOError(std::string("s3://") + "user:password@" + "bucket/key?signature=private-value");
      }
      return reader_.Read(std::min<int64_t>(nbytes, 1), out);
    }

    arrow::Result<std::shared_ptr<arrow::Buffer>> Read(int64_t nbytes) override
    {
      return reader_.Read(std::min<int64_t>(nbytes, 1));
    }

    arrow::Status Close() override { return reader_.Close(); }
    bool closed() const override { return reader_.closed(); }
    arrow::Result<int64_t> Tell() const override { return reader_.Tell(); }

  private:
    arrow::io::BufferReader reader_;
    int64_t fail_after_;
  };

  std::string readAll(ArrowInputStream& stream)
  {
    std::array<XMLByte, 7> buffer{};
    std::string output;
    XMLSize_t count;
    while ((count = stream.readBytes(buffer.data(), buffer.size())) != 0)
    {
      output.append(reinterpret_cast<const char*>(buffer.data()), count);
    }
    return output;
  }
}

START_TEST(S3InputSource, "$Id$")

START_SECTION(([EXTRA] partial reads, zero reads, EOF, ownership and positions))
{
  for (const std::string data : {"", "x", "BZ", "<root>remote XML</root>"})
  {
    auto input = std::make_shared<ChunkedStream>(data);
    {
      ArrowInputStream stream(input);
      TEST_EQUAL(stream.curPos(), 0)
      TEST_EQUAL(stream.readBytes(nullptr, 0), 0)
      TEST_EQUAL(readAll(stream), data)
      TEST_EQUAL(stream.curPos(), data.size())
      XMLByte byte;
      TEST_EQUAL(stream.readBytes(&byte, 1), 0)
      TEST_EQUAL(stream.readBytes(&byte, 1), 0)
      TEST_TRUE(stream.getContentType() == nullptr)
    }
    TEST_TRUE(input->closed())
  }
}
END_SECTION

START_SECTION(([EXTRA] bounded reads and null arguments))
{
  auto input = std::make_shared<arrow::io::BufferReader>(arrow::Buffer::FromString(std::string(200000, 'x')));
  ArrowInputStream stream(input);
  std::array<XMLByte, 100000> buffer{};
  TEST_EQUAL(stream.readBytes(buffer.data(), buffer.size()), 3)
  TEST_EQUAL(stream.readBytes(buffer.data(), buffer.size()), 65536)
  TEST_EQUAL(stream.curPos(), 65539)
  TEST_EXCEPTION(Exception::IOException, stream.readBytes(nullptr, 1))
  TEST_EXCEPTION(Exception::IOException, (void)ArrowInputStream(nullptr))
}
END_SECTION

START_SECTION(([EXTRA] read errors are exceptions, not EOF))
{
  auto input = std::make_shared<ChunkedStream>("abcdef", 3);
  ArrowInputStream stream(input);
  XMLByte buffer[8];
  TEST_EQUAL(stream.readBytes(buffer, sizeof(buffer)), 3)
  TEST_EXCEPTION(Exception::IOException, stream.readBytes(buffer, sizeof(buffer)))
  TEST_EQUAL(stream.curPos(), 3)
  try
  {
    stream.readBytes(buffer, sizeof(buffer));
    TEST_TRUE(false)
  }
  catch (const Exception::IOException& error)
  {
    const std::string message = error.getMessage();
    TEST_TRUE(message.find("read/decompression") != std::string::npos)
    TEST_TRUE(message.find("IOError") != std::string::npos)
    TEST_TRUE(message.find("password") == std::string::npos)
    TEST_TRUE(message.find("signature") == std::string::npos)
    TEST_TRUE(message.find("private-value") == std::string::npos)
  }
  auto failing = std::make_shared<ChunkedStream>("abcdef", 1);
  TEST_EXCEPTION(Exception::IOException, (void)ArrowInputStream(failing))
  TEST_TRUE(failing->closed())
}
END_SECTION

START_SECTION(([EXTRA] gzip and bzip2 magic detection, EOF and truncation))
{
  const std::string gzip(
    "\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\x03\xb3\x29\xca\xcf\x2f\xb1\x2b\x4a\xcd\xcd"
    "\x2f\x49\x55\x88\xf0\xf5\xb1\xd1\x07\x0b\x00\x00\x53\x37\x5f\xa3\x17\x00\x00\x00", 40);
  const std::string bzip2(
    "\x42\x5a\x68\x39\x31\x41\x59\x26\x53\x59\xc6\x41\xf5\xea\x00\x00\x01\x9f\x80\x40"
    "\x00\x80\x05\x00\x06\x00\x40\x02\x02\x94\x00\x20\x00\x31\x03\x40\xd0\x1a\x86\x9a"
    "\x79\x47\xea\x8e\xae\xc4\xe3\xa2\x65\x08\x81\xbc\x7c\x5d\xc9\x14\xe1\x42\x43\x19"
    "\x07\xd7\xa8", 63);
  for (const auto& compressed : {gzip, bzip2})
  {
    ArrowInputStream stream(std::make_shared<ChunkedStream>(compressed));
    TEST_EQUAL(readAll(stream), "<root>remote XML</root>")
    TEST_EQUAL(stream.curPos(), 23)
    XMLByte byte;
    TEST_EQUAL(stream.readBytes(&byte, 1), 0)
    ArrowInputStream truncated(std::make_shared<ChunkedStream>(compressed.substr(0, compressed.size() - 5)));
    TEST_EXCEPTION(Exception::IOException, readAll(truncated))
    auto corrupt = compressed;
    corrupt[10] ^= 0x7f;
    ArrowInputStream damaged(std::make_shared<ChunkedStream>(corrupt));
    TEST_EXCEPTION(Exception::IOException, readAll(damaged))
  }
}
END_SECTION

START_SECTION(([EXTRA] input source construction does not access S3))
{
  Internal::xercesInitialize();
  S3InputSource source(std::string("s3://") + "user:password@" + "unused-bucket/example.mzML?signature=private-value#private-fragment");
  TEST_TRUE(source.getSystemId() != nullptr)
  TEST_TRUE(xercesc::XMLString::equals(source.getSystemId(), u"s3://unused-bucket/example.mzML"))
  TEST_EXCEPTION(Exception::ParseError, (void)S3InputSource("file:///unused.xml"))
  S3InputSource malformed(std::string("s3://") + "user:password@" + "[invalid/key?signature=private-value");
  try
  {
    std::unique_ptr<xercesc::BinInputStream> stream(malformed.makeStream());
    TEST_TRUE(false)
  }
  catch (const Exception::IOException& error)
  {
    const std::string message = error.getMessage();
    TEST_TRUE(message.find("URI resolution") != std::string::npos)
    TEST_TRUE(message.find("password") == std::string::npos)
    TEST_TRUE(message.find("signature") == std::string::npos)
    TEST_TRUE(message.find("private-value") == std::string::npos)
  }
}
END_SECTION

END_TEST
