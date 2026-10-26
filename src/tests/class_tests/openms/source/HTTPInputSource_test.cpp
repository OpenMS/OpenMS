// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h>
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
  // A file:// URL for a local path; libcurl reads it like a remote resource, so no network is needed
  std::string fileURL(std::string path)
  {
    std::replace(path.begin(), path.end(), '\\', '/');
    std::string url = path.empty() || path[0] != '/' ? "file:///" : "file://"; // Windows: file:///C:/...
    for (const char c : path)
    {
      if (c == ' ') { url += "%20"; }
      else if (c == '%') { url += "%25"; }
      else if (c == '#') { url += "%23"; }
      else { url += c; }
    }
    return url;
  }

  std::string readFile(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }

  std::string readAll(InputStream& stream)
  {
    std::string result;
    unsigned char buffer[1000];
    while (const Size count = stream.readBytes(buffer, sizeof(buffer)))
    {
      result.append(reinterpret_cast<const char*>(buffer), count);
    }
    return result;
  }
}

START_TEST(HTTPInputSource, "$Id$")

START_SECTION((static bool isURL(const std::string& url)))
{
  TEST_TRUE(HTTPInputSource::isURL("http://example.org/run.mzML"))
  TEST_TRUE(HTTPInputSource::isURL("https://example.org/run.mzML.gz?token=1"))
  TEST_TRUE(HTTPInputSource::isURL("ftp://ftp.example.org/pride/run.mzML"))
  TEST_FALSE(HTTPInputSource::isURL("s3://bucket/run.mzML"))
  TEST_FALSE(HTTPInputSource::isURL("file:///data/run.mzML")) // local files are read directly
  TEST_FALSE(HTTPInputSource::isURL("data/http://run.mzML"))
  TEST_FALSE(HTTPInputSource::isURL("C:/data/run.mzML"))
}
END_SECTION

START_SECTION((static std::string redactURL(const std::string& url)))
{
  // user info, query and fragment are removed
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("https://user:secret@example.org:8443/data/run.mzML?X-Amz-Signature=abc#top"),
                    "https://example.org:8443/data/run.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("ftp://anonymous@ftp.example.org/pub/run.mzML"), "ftp://ftp.example.org/pub/run.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("http://example.org?token=1"), "http://example.org")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("s3://key:secret@bucket/run.mzML"), "s3://bucket/run.mzML")
  // an '@' in the path is no user info
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("http://example.org/a@b/run.mzML"), "http://example.org/a@b/run.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("file:///data/run.mzML"), "file:///data/run.mzML")
  // local file names are unchanged
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("/data/run?1#2@x.mzML"), "/data/run?1#2@x.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("C:\\data\\run#1.mzML"), "C:\\data\\run#1.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("C://data/run?1.mzML"), "C://data/run?1.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL("data/http://run?1.mzML"), "data/http://run?1.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURL(""), "")
}
END_SECTION

START_SECTION((explicit HTTPInputSource(const std::string& url)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, HTTPInputSource("s3://bucket/run.mzML"))
  TEST_EXCEPTION(Exception::IllegalArgument, HTTPInputSource(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML")))
  HTTPInputSource source("http://example.org/run.mzML"); // nothing is transferred yet
  NOT_TESTABLE
}
END_SECTION

START_SECTION((std::unique_ptr<InputStream> makeStream() const))
{
  // plain content
  const std::string plain_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML");
  auto plain = HTTPInputSource(fileURL(plain_path)).makeStream();
  TEST_EQUAL(plain->curPos(), 0)
  const std::string plain_content = readAll(*plain);
  TEST_EQUAL(plain_content.empty(), false)
  TEST_TRUE(plain_content == readFile(plain_path))
  TEST_EQUAL(plain->curPos(), plain_content.size())

  // gzip and bzip2 compressed content is decompressed
  const std::string gz_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz");
  GzipInputStream gzip(gz_path);
  auto gz = HTTPInputSource(fileURL(gz_path)).makeStream();
  TEST_TRUE(readAll(*gz) == readAll(gzip))
  const std::string bz2_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2");
  Bzip2InputStream bzip2(bz2_path);
  auto bz2 = HTTPInputSource(fileURL(bz2_path)).makeStream();
  TEST_TRUE(readAll(*bz2) == readAll(bzip2))

  // a missing resource, a server that cannot be reached (port 1 on the local host)
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource(fileURL(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist.mzML"))).makeStream())
  TEST_EXCEPTION(Exception::IOException, HTTPInputSource("http://127.0.0.1:1/run.mzML").makeStream())
  // errors do not name the credentials of the URL
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IOException, HTTPInputSource("http://user:secret@127.0.0.1:1/run.mzML?token=1").makeStream(),
                              "IO error for file 'http://127.0.0.1:1/run.mzML'")
}
END_SECTION

END_TEST
