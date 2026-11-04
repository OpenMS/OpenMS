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

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <thread>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX // std::min, not the min macro of windows.h
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "Ws2_32.lib")
#else
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>
#endif

using namespace OpenMS;

namespace
{
#ifdef _WIN32
  using Socket = SOCKET;
  const Socket NO_SOCKET = INVALID_SOCKET;
  void closeSocket(Socket s) { closesocket(s); }
#else
  using Socket = int;
  const Socket NO_SOCKET = -1;
  void closeSocket(Socket s) { close(s); }
#endif

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

  /**
    A minimal HTTP/1.0 server on 127.0.0.1 for the tests: GET and HEAD of the given files, byte ranges (206; 416 beyond
    the end), 404 for other paths. Paths under /norange/ serve the same files but ignore range requests (200 with the
    whole file), as some servers do. Paths under /encoded/ serve the gzip file of the path (path + ".gz") with gzip
    transfer compression (Content-Encoding), if the client accepts it (otherwise 406). Paths under /moved/ redirect to
    the path without /moved/. Paths under /signed/ act like a URL presigned for GET: only a GET with a query (the
    signature) is served, anything else is 403.
  */
  class TestServer
  {
public:
    explicit TestServer(std::map<std::string, std::string> files) :
      files_(std::move(files))
    {
#ifdef _WIN32
      WSADATA data;
      WSAStartup(MAKEWORD(2, 2), &data);
#endif
      listener_ = socket(AF_INET, SOCK_STREAM, 0);
      sockaddr_in address{};
      address.sin_family = AF_INET;
      address.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
      address.sin_port = 0; // any free port
      bind(listener_, reinterpret_cast<sockaddr*>(&address), sizeof(address));
      listen(listener_, 16);
      socklen_t length = sizeof(address);
      getsockname(listener_, reinterpret_cast<sockaddr*>(&address), &length);
      port_ = ntohs(address.sin_port);
      worker_ = std::thread([this] { run_(); });
    }

    ~TestServer()
    {
      stop_ = true;
      // wake up accept()
      Socket wake = socket(AF_INET, SOCK_STREAM, 0);
      sockaddr_in address{};
      address.sin_family = AF_INET;
      address.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
      address.sin_port = htons(port_);
      connect(wake, reinterpret_cast<sockaddr*>(&address), sizeof(address));
      closeSocket(wake);
      worker_.join();
      closeSocket(listener_);
#ifdef _WIN32
      WSACleanup();
#endif
    }

    std::string url(const std::string& path) const
    {
      return "http://127.0.0.1:" + std::to_string(port_) + path;
    }

    unsigned short port() const { return port_; }

private:
    void run_()
    {
      while (true)
      {
        Socket connection = accept(listener_, nullptr, nullptr);
        if (stop_)
        {
          if (connection != NO_SOCKET) { closeSocket(connection); }
          return;
        }
        if (connection == NO_SOCKET) { continue; }
        serve_(connection);
        closeSocket(connection);
      }
    }

    void serve_(Socket connection)
    {
      std::string request;
      char buffer[4096];
      while (request.find("\r\n\r\n") == std::string::npos)
      {
        const auto count = recv(connection, buffer, sizeof(buffer), 0);
        if (count <= 0) { return; }
        request.append(buffer, static_cast<size_t>(count));
      }
      std::istringstream lines(request);
      std::string method, path;
      lines >> method >> path;
      const bool query = path.find('?') != std::string::npos;
      path = path.substr(0, path.find('?')); // the query does not select the file
      std::string lowered = request;
      std::transform(lowered.begin(), lowered.end(), lowered.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
      if (path.rfind("/moved/", 0) == 0)
      {
        // with a body, which is not the resource
        const std::string body = "moved";
        send_(connection, "HTTP/1.0 302 Found\r\nLocation: " + path.substr(6) + "\r\nContent-Length: " + std::to_string(body.size()) + "\r\n\r\n"
                          + (method == "HEAD" ? "" : body));
        return;
      }
      if (path.rfind("/signed/", 0) == 0)
      {
        if (method != "GET" || !query)
        {
          send_(connection, "HTTP/1.0 403 Forbidden\r\nContent-Length: 0\r\n\r\n");
          return;
        }
        path = path.substr(7);
      }
      bool ranges = true;
      if (path.rfind("/norange/", 0) == 0)
      {
        path = path.substr(8);
        ranges = false;
      }
      std::string headers;
      if (path.rfind("/encoded/", 0) == 0)
      {
        const auto accept = lowered.find("\r\naccept-encoding:");
        if (accept == std::string::npos || lowered.substr(accept + 2, lowered.find("\r\n", accept + 2) - accept - 2).find("gzip") == std::string::npos)
        {
          send_(connection, "HTTP/1.0 406 Not Acceptable\r\nContent-Length: 0\r\n\r\n");
          return;
        }
        path = path.substr(8) + ".gz";
        headers = "Content-Encoding: gzip\r\n";
        ranges = false;
      }
      // the requested byte range, if any ("Range: bytes=first-last")
      long long first = -1, last = -1;
      const auto range = lowered.find("\r\nrange: bytes=");
      if (ranges && range != std::string::npos)
      {
        std::istringstream spec(lowered.substr(range + 15));
        char dash = 0;
        spec >> first >> dash >> last;
      }
      const auto file = files_.find(path);
      if (file == files_.end())
      {
        send_(connection, "HTTP/1.0 404 Not Found\r\nContent-Length: 0\r\n\r\n");
        return;
      }
      const std::string& content = file->second;
      const long long size = static_cast<long long>(content.size());
      std::string status = "200 OK", body = content;
      if (first >= 0)
      {
        if (first >= size)
        {
          send_(connection, "HTTP/1.0 416 Range Not Satisfiable\r\nContent-Range: bytes */" + std::to_string(size) + "\r\nContent-Length: 0\r\n\r\n");
          return;
        }
        last = std::min(last, size - 1);
        status = "206 Partial Content";
        body = content.substr(static_cast<size_t>(first), static_cast<size_t>(last - first + 1));
        headers += "Content-Range: bytes " + std::to_string(first) + "-" + std::to_string(last) + "/" + std::to_string(size) + "\r\n";
      }
      std::string response = "HTTP/1.0 " + status + "\r\nContent-Length: " + std::to_string(body.size()) + "\r\n" + headers + "\r\n";
      if (method != "HEAD") { response += body; }
      send_(connection, response);
    }

    static void send_(Socket connection, const std::string& data)
    {
      int flags = 0;
#ifdef MSG_NOSIGNAL
      flags = MSG_NOSIGNAL; // the client may close early (e.g. after the range it needs)
#endif
      size_t sent = 0;
      while (sent < data.size())
      {
        const auto count = ::send(connection, data.data() + sent, static_cast<int>(data.size() - sent), flags);
        if (count <= 0) { return; }
        sent += static_cast<size_t>(count);
      }
    }

    std::map<std::string, std::string> files_;
    Socket listener_ = NO_SOCKET;
    unsigned short port_ = 0;
    std::atomic<bool> stop_{false};
    std::thread worker_;
  };
}

START_TEST(HTTPInputSource, "$Id$")

#ifndef _WIN32
// macOS has no MSG_NOSIGNAL: a client closing early must not kill the test
signal(SIGPIPE, SIG_IGN);
#endif

// requests to the test server must not go through a proxy of the environment
#ifdef _WIN32
_putenv_s("no_proxy", "127.0.0.1");
_putenv_s("NO_PROXY", "127.0.0.1");
#else
setenv("no_proxy", "127.0.0.1", 1);
setenv("NO_PROXY", "127.0.0.1", 1);
#endif

const std::string served_plain_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML");
const std::string served_indexed_path = OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML");
const std::string served_plain = readFile(served_plain_path);
// the files the test server serves at these paths
const std::map<std::string, std::string> local_paths{{"/run.mzML", served_plain_path},
                                                     {"/run.mzML.gz", OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz")},
                                                     {"/run.mzML.bz2", OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2")},
                                                     // gzip of the compressed files, for /encoded/run.mzML.gz and /encoded/run.mzML.bz2
                                                     {"/run.mzML.gz.gz", OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz.gz")},
                                                     {"/run.mzML.bz2.gz", OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2.gz")},
                                                     {"/indexed.mzML", served_indexed_path}};
std::map<std::string, std::string> served;
for (const auto& [path, local] : local_paths) { served[path] = readFile(local); }
// a binary array declared zlib-compressed that is not: the handler throws while decoding it
std::string corrupt = served["/run.mzML"];
const std::string no_compression = "accession=\"MS:1000576\" name=\"no compression\"";
corrupt.replace(corrupt.find(no_compression), no_compression.size(), "accession=\"MS:1000574\" name=\"zlib compression\"");
served["/corrupt.mzML"] = corrupt;
served["/empty.mzML"] = "";
// a gzip header, then no valid deflate data (for /encoded/broken.mzML)
served["/broken.mzML.gz"] = std::string("\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03", 10) + std::string(100, '\xff');
TestServer server(served);

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

START_SECTION((static std::string redactURLs(const std::string& text)))
{
  // e.g. Arrow's diagnostic for an invalid port, which repeats the whole URI
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("Invalid: Cannot parse URI: 's3://KEY:SECRET@bucket:999999/run.mzML?region=us-east-1&token=TOKEN'"),
                    "Invalid: Cannot parse URI: 's3://bucket:999999/run.mzML'")
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("from ftp://user:pw@host/a.mzML to https://example.org/b?sig=1 failed"),
                    "from ftp://host/a.mzML to https://example.org/b failed")
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("no URL here"), "no URL here")
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs(""), "")
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("://nothing before"), "://nothing before")
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("file C://data/run?1.mzML"), "file C://data/run?1.mzML")
}
END_SECTION

START_SECTION((static std::string redactURLs(const std::string& text, const std::string& url)))
{
  const std::string url = "s3://KEY:SE%2FCRET@bucket/run.mzML?region=eu-1&scheme=MY_SCHEME#part";
  // Arrow repeats a query value on its own
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("Invalid S3 connection scheme 'MY_SCHEME'", url), "Invalid S3 connection scheme '***'")
  // user name and password (as given and decoded), query values, fragment; URLs as by redactURLs(text)
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("KEY, SE/CRET, SE%2FCRET, eu-1, part: s3://KEY:SE%2FCRET@bucket/run.mzML?scheme=MY_SCHEME", url),
                    "***, ***, ***, ***, ***: s3://bucket/run.mzML")
  // only whole words
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("port 1234, value 1, KEYS", "http://KEY@example.org/x?a=1"), "port 1234, value ***, KEYS")
  // nothing to remove: no such parts, a local file name
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("region eu-1 of s3://bucket/run.mzML", "s3://bucket/run.mzML"), "region eu-1 of s3://bucket/run.mzML")
  TEST_STRING_EQUAL(HTTPInputSource::redactURLs("cannot read /data/run?1#2.mzML", "/data/run?1#2.mzML"), "cannot read /data/run?1#2.mzML")
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
  TEST_FALSE(plain_content.empty())
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

START_SECTION(([EXTRA] makeStream from an HTTP server))
{
  auto plain = HTTPInputSource(server.url("/run.mzML")).makeStream();
  TEST_TRUE(readAll(*plain) == served_plain)
  // gzip and bzip2 compressed content is decompressed as by the local readers
  GzipInputStream gzip(local_paths.at("/run.mzML.gz"));
  auto gz = HTTPInputSource(server.url("/run.mzML.gz")).makeStream();
  TEST_TRUE(readAll(*gz) == readAll(gzip))
  Bzip2InputStream bzip2(local_paths.at("/run.mzML.bz2"));
  auto bz2 = HTTPInputSource(server.url("/run.mzML.bz2")).makeStream();
  TEST_TRUE(readAll(*bz2) == readAll(bzip2))
  // gzip transfer compression is requested
  GzipInputStream encoded_gzip(local_paths.at("/run.mzML.gz"));
  auto encoded = HTTPInputSource(server.url("/encoded/run.mzML")).makeStream();
  TEST_TRUE(readAll(*encoded) == readAll(encoded_gzip))
  // ... and decompressed by OpenMS as the stream is read, not by libcurl (which decodes all it receives at once, too
  // much for a highly compressed response): corrupt data is a ParseError, as for a compressed resource
  TEST_EXCEPTION(Exception::ParseError, readAll(*HTTPInputSource(server.url("/encoded/broken.mzML")).makeStream()))
  // ... around a compressed resource: both are decompressed
  GzipInputStream encoded_gz_content(local_paths.at("/run.mzML.gz"));
  auto encoded_gz = HTTPInputSource(server.url("/encoded/run.mzML.gz")).makeStream();
  TEST_TRUE(readAll(*encoded_gz) == readAll(encoded_gz_content))
  Bzip2InputStream encoded_bz2_content(local_paths.at("/run.mzML.bz2"));
  auto encoded_bz2 = HTTPInputSource(server.url("/encoded/run.mzML.bz2")).makeStream();
  TEST_TRUE(readAll(*encoded_bz2) == readAll(encoded_bz2_content))
  // a redirect is followed (its body is not the content)
  auto moved = HTTPInputSource(server.url("/moved/run.mzML")).makeStream();
  TEST_TRUE(readAll(*moved) == served_plain)
  // HTTP 404
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource(server.url("/missing.mzML")).makeStream())
}
END_SECTION

START_SECTION(([EXTRA] MzMLFile::load from an HTTP server))
{
  for (const std::string path : {"/run.mzML", "/run.mzML.gz", "/run.mzML.bz2", "/indexed.mzML", "/encoded/run.mzML", "/encoded/run.mzML.gz", "/encoded/run.mzML.bz2"})
  {
    PeakMap local, remote;
    // gzip transfer compression of /encoded/...: of /run.mzML.gz for /encoded/run.mzML, of the compressed file otherwise
    const bool encoded = path.rfind("/encoded/", 0) == 0;
    MzMLFile().load(local_paths.at(path == "/encoded/run.mzML" ? "/run.mzML.gz" : (encoded ? path.substr(8) : path)), local);
    MzMLFile().load(server.url(path), remote);
    TEST_EQUAL(remote.size(), local.size())
    TEST_TRUE(remote.getSpectra() == local.getSpectra())
    TEST_TRUE(remote.getChromatograms() == local.getChromatograms())
    TEST_STRING_EQUAL(remote.getLoadedFilePath(), server.url(path)) // the URL as given
  }
  PeakMap missing;
  TEST_EXCEPTION(Exception::ParseError, MzMLFile().load(server.url("/missing.mzML"), missing))
}
END_SECTION

START_SECTION((static Int64 objectSize(const std::string& url)))
{
  const Int64 size = static_cast<Int64>(served_plain.size());
  TEST_EQUAL(HTTPInputSource::objectSize(server.url("/run.mzML")), size)
  // a GET of the first byte, not a HEAD request, which a URL presigned for GET refuses
  TEST_EQUAL(HTTPInputSource::objectSize(server.url("/signed/run.mzML?X-Amz-Signature=abc")), size)
  // a server that ignores the range (the whole resource is not read), a redirect, an empty resource (416)
  TEST_EQUAL(HTTPInputSource::objectSize(server.url("/norange/run.mzML")), size)
  TEST_EQUAL(HTTPInputSource::objectSize(server.url("/moved/run.mzML")), size)
  TEST_EQUAL(HTTPInputSource::objectSize(server.url("/empty.mzML")), 0)
  TEST_EQUAL(HTTPInputSource::objectSize(fileURL(served_plain_path)), size)
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource::objectSize(server.url("/missing.mzML")))
  TEST_EXCEPTION(Exception::IOException, HTTPInputSource::objectSize(server.url("/signed/run.mzML"))) // HTTP 403
  TEST_EXCEPTION(Exception::IllegalArgument, HTTPInputSource::objectSize(served_plain_path))
}
END_SECTION

START_SECTION((static std::string readRange(const std::string& url, Int64 offset, Int64 length)))
{
  const Int64 size = static_cast<Int64>(served_plain.size());
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), 100, 50) == served_plain.substr(100, 50))
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), 0, 10) == served_plain.substr(0, 10))
  TEST_TRUE(HTTPInputSource::readRange(fileURL(served_plain_path), 100, 50) == served_plain.substr(100, 50))
  // fewer bytes at the end
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), size - 5, 100) == served_plain.substr(served_plain.size() - 5))
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), 100, 0).empty())
  TEST_TRUE(HTTPInputSource::readRange(server.url("/signed/run.mzML?X-Amz-Signature=abc"), 100, 50) == served_plain.substr(100, 50))
  // a server that ignores the range
  TEST_TRUE(HTTPInputSource::readRange(server.url("/norange/run.mzML"), 100, 50) == served_plain.substr(100, 50))
  TEST_TRUE(HTTPInputSource::readRange(server.url("/norange/run.mzML"), size - 5, 100) == served_plain.substr(served_plain.size() - 5))
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource::readRange(server.url("/missing.mzML"), 0, 10))
}
END_SECTION

START_SECTION(([EXTRA] IndexedMzMLDecoder reads the index over HTTP))
{
  IndexedMzMLDecoder decoder;
  const std::streampos local_offset = decoder.findIndexListOffset(served_indexed_path);
  const std::streampos remote_offset = decoder.findIndexListOffset(server.url("/indexed.mzML"));
  TEST_EQUAL(remote_offset, local_offset)
  IndexedMzMLDecoder::OffsetVector local_spectra, local_chromatograms, remote_spectra, remote_chromatograms;
  TEST_EQUAL(decoder.parseOffsets(served_indexed_path, local_offset, local_spectra, local_chromatograms), 0)
  TEST_EQUAL(decoder.parseOffsets(server.url("/indexed.mzML"), remote_offset, remote_spectra, remote_chromatograms), 0)
  TEST_TRUE(remote_spectra == local_spectra)
  TEST_TRUE(remote_chromatograms == local_chromatograms)
}
END_SECTION

START_SECTION(([EXTRA] errors while decoding remote binary data do not show credentials))
{
  const std::string url = "http://user:secret@127.0.0.1:" + std::to_string(server.port()) + "/corrupt.mzML?token=TOKEN";
  PeakMap map;
  std::string message;
  try
  {
    MzMLFile().load(url, map);
  }
  catch (const Exception::ParseError& e)
  {
    message = e.what();
  }
  TEST_EQUAL(message.empty(), false)
  TEST_EQUAL(message.find("secret"), std::string::npos)
  TEST_EQUAL(message.find("TOKEN"), std::string::npos)
}
END_SECTION

END_TEST
