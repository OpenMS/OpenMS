// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <csignal>
#include <cstdlib>
#include <fstream>
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

  std::string readFile(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream content;
    content << in.rdbuf();
    return content.str();
  }

  std::string readAll(InputStream& stream)
  {
    std::string content;
    unsigned char buffer[1000];
    while (Size count = stream.readBytes(buffer, sizeof(buffer)))
    {
      content.append(reinterpret_cast<const char*>(buffer), count);
    }
    return content;
  }

  /**
    A minimal HTTP/1.0 server on 127.0.0.1 for the tests: GET and HEAD of the given files, byte ranges (206), 404 for
    other paths. Paths under /norange/ serve the same files but ignore range requests (200 with the whole file), as
    some servers do.
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
      bool ranges = true;
      if (path.rfind("/norange/", 0) == 0)
      {
        path = path.substr(8);
        ranges = false;
      }
      // the requested byte range, if any ("Range: bytes=first-last")
      long long first = -1, last = -1;
      std::string lowered = request;
      std::transform(lowered.begin(), lowered.end(), lowered.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
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
      std::string status = "200 OK", body = content, headers;
      if (first >= 0)
      {
        if (first >= size)
        {
          send_(connection, "HTTP/1.0 416 Range Not Satisfiable\r\nContent-Length: 0\r\n\r\n");
          return;
        }
        last = std::min(last, size - 1);
        status = "206 Partial Content";
        body = content.substr(static_cast<size_t>(first), static_cast<size_t>(last - first + 1));
        headers = "Content-Range: bytes " + std::to_string(first) + "-" + std::to_string(last) + "/" + std::to_string(size) + "\r\n";
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

const std::string plain_path = OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML");
const std::string indexed_path = OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML");
const std::string plain = readFile(plain_path);
const std::string indexed = readFile(indexed_path);
// the files the server serves at these paths
const std::map<std::string, std::string> local_paths{{"/run.mzML", plain_path},
                                                     {"/run.mzML.gz", OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.gz")},
                                                     {"/run.mzML.bz2", OPENMS_GET_TEST_DATA_PATH("MzMLFile_6_uncompressed.mzML.bz2")}};
TestServer server({{"/run.mzML", plain},
                   {"/run.mzML.gz", readFile(local_paths.at("/run.mzML.gz"))},
                   {"/run.mzML.bz2", readFile(local_paths.at("/run.mzML.bz2"))},
                   {"/indexed.mzML", indexed}});

START_SECTION((static bool isHTTPURI(const std::string& url)))
{
  TEST_TRUE(HTTPInputSource::isHTTPURI("http://example.org/run.mzML"))
  TEST_TRUE(HTTPInputSource::isHTTPURI("https://example.org/run.mzML?token=1"))
  TEST_TRUE(HTTPInputSource::isHTTPURI("ftp://ftp.pride.ebi.ac.uk/run.mzML"))
  TEST_TRUE(HTTPInputSource::isHTTPURI("HTTPS://example.org/run.mzML"))
  TEST_FALSE(HTTPInputSource::isHTTPURI("s3://bucket/run.mzML"))
  TEST_FALSE(HTTPInputSource::isHTTPURI("data/http://run.mzML"))
  TEST_FALSE(HTTPInputSource::isHTTPURI(plain_path))
}
END_SECTION

START_SECTION((explicit HTTPInputSource(const std::string& url)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, HTTPInputSource("s3://bucket/run.mzML"))
  TEST_EXCEPTION(Exception::IllegalArgument, HTTPInputSource{plain_path})
  TEST_EXCEPTION(Exception::IllegalArgument, HTTPInputSource::objectSize(plain_path))
}
END_SECTION

START_SECTION((std::unique_ptr<InputStream> makeStream() const))
{
  // the stream is the file's content, decompressed as for a local file
  {
    auto stream = HTTPInputSource(server.url("/run.mzML")).makeStream();
    TEST_TRUE(readAll(*stream) == plain)
    TEST_EQUAL(stream->curPos(), plain.size())
  }
  for (const std::string path : {"/run.mzML.gz", "/run.mzML.bz2"})
  {
    const std::string& local_path = local_paths.at(path);
    std::unique_ptr<InputStream> local(CompressedInputSource(local_path, path.substr(path.size() - 3) == "bz2" ? "BZ" : "\x1f\x8b").makeStream());
    const std::string expected = readAll(*local);
    auto stream = HTTPInputSource(server.url(path)).makeStream();
    TEST_TRUE(readAll(*stream) == expected)
    TEST_EQUAL(stream->curPos(), expected.size())
  }
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource(server.url("/missing.mzML")).makeStream())
}
END_SECTION

START_SECTION(([EXTRA] MzMLFile::load from HTTP URLs))
{
  for (const std::string path : {"/run.mzML", "/run.mzML.gz", "/run.mzML.bz2"})
  {
    PeakMap local, remote;
    MzMLFile().load(local_paths.at(path), local);
    MzMLFile().load(server.url(path), remote);
    TEST_EQUAL(remote.size(), local.size())
    TEST_TRUE(remote.getSpectra() == local.getSpectra())
    TEST_TRUE(remote.getChromatograms() == local.getChromatograms())
    // the URL is kept as it is
    TEST_STRING_EQUAL(remote.getLoadedFilePath(), server.url(path))
  }
  // an indexed file, whose index is not used for a remote load
  PeakMap local_indexed, remote_indexed;
  MzMLFile().load(indexed_path, local_indexed);
  MzMLFile().load(server.url("/indexed.mzML"), remote_indexed);
  TEST_TRUE(remote_indexed.getSpectra() == local_indexed.getSpectra())
  PeakMap missing;
  TEST_EXCEPTION(Exception::ParseError, MzMLFile().load(server.url("/missing.mzML"), missing))
}
END_SECTION

START_SECTION((static Int64 objectSize(const std::string& url)))
{
  TEST_EQUAL(HTTPInputSource::objectSize(server.url("/run.mzML")), static_cast<Int64>(plain.size()))
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource::objectSize(server.url("/missing.mzML")))
}
END_SECTION

START_SECTION((static std::string readRange(const std::string& url, Int64 offset, Int64 length)))
{
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), 100, 50) == plain.substr(100, 50))
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), 0, 10) == plain.substr(0, 10))
  // fewer bytes at the end
  const Int64 size = static_cast<Int64>(plain.size());
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), size - 5, 100) == plain.substr(plain.size() - 5))
  TEST_TRUE(HTTPInputSource::readRange(server.url("/run.mzML"), 100, 0).empty())
  // a server that ignores the range
  TEST_TRUE(HTTPInputSource::readRange(server.url("/norange/run.mzML"), 100, 50) == plain.substr(100, 50))
  TEST_TRUE(HTTPInputSource::readRange(server.url("/norange/run.mzML"), size - 5, 100) == plain.substr(plain.size() - 5))
  TEST_EXCEPTION(Exception::FileNotFound, HTTPInputSource::readRange(server.url("/missing.mzML"), 0, 10))
}
END_SECTION

START_SECTION(([EXTRA] IndexedMzMLDecoder reads the index of an HTTP URL))
{
  IndexedMzMLDecoder decoder;
  const std::streampos local_offset = decoder.findIndexListOffset(indexed_path);
  const std::streampos remote_offset = decoder.findIndexListOffset(server.url("/indexed.mzML"));
  TEST_EQUAL(remote_offset, local_offset)
  IndexedMzMLDecoder::OffsetVector local_spectra, local_chromatograms, remote_spectra, remote_chromatograms;
  TEST_EQUAL(decoder.parseOffsets(indexed_path, local_offset, local_spectra, local_chromatograms), 0)
  TEST_EQUAL(decoder.parseOffsets(server.url("/indexed.mzML"), remote_offset, remote_spectra, remote_chromatograms), 0)
  TEST_TRUE(remote_spectra == local_spectra)
  TEST_TRUE(remote_chromatograms == local_chromatograms)
}
END_SECTION

START_SECTION(([EXTRA] a server that cannot be reached))
{
  unsigned short closed_port = 0;
  {
    TestServer stopped({});
    closed_port = stopped.port();
  }
  const std::string url = "http://127.0.0.1:" + std::to_string(closed_port) + "/run.mzML";
  TEST_EXCEPTION(Exception::IOException, HTTPInputSource(url).makeStream())
}
END_SECTION

END_TEST
