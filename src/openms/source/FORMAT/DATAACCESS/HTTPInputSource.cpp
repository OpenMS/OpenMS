// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/DecompressingInputStream.h>
#include <OpenMS/SYSTEM/CurlInit.h>

#include <curl/curl.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <optional>
#include <string_view>
#include <vector>

namespace OpenMS
{
  namespace
  {
#if LIBCURL_VERSION_NUM >= 0x075500 // 7.85.0
    // The protocols of @p wanted that this libcurl supports, comma-separated (CURLOPT_PROTOCOLS_STR rejects others)
    std::string supportedProtocols(std::initializer_list<const char*> wanted)
    {
      const curl_version_info_data* info = curl_version_info(CURLVERSION_NOW);
      std::string list;
      for (const char* protocol : wanted)
      {
        for (const char* const* supported = info->protocols; supported != nullptr && *supported != nullptr; ++supported)
        {
          if (std::strcmp(*supported, protocol) == 0)
          {
            list += (list.empty() ? "" : ",") + std::string(protocol);
            break;
          }
        }
      }
      return list;
    }
#endif

    // @p text with its %XX escapes decoded
    std::string percentDecoded(std::string_view text)
    {
      std::string result;
      for (Size i = 0; i < text.size(); ++i)
      {
        if (text[i] == '%' && i + 2 < text.size() && std::isxdigit(static_cast<unsigned char>(text[i + 1])) && std::isxdigit(static_cast<unsigned char>(text[i + 2])))
        {
          result += static_cast<char>(std::stoi(std::string(text.substr(i + 1, 2)), nullptr, 16));
          i += 2;
        }
        else
        {
          result += text[i];
        }
      }
      return result;
    }

    // The options all transfers share; false if the URL or the protocol restriction cannot be set (the transfer must not
    // start then). @p error receives libcurl's error messages.
    bool setUpTransfer(CURL* easy, const std::string& url, char* error)
    {
      static const std::string user_agent = "OpenMS/" + VersionInfo::getVersion();
      bool ok = curl_easy_setopt(easy, CURLOPT_URL, url.c_str()) == CURLE_OK;
      curl_easy_setopt(easy, CURLOPT_ERRORBUFFER, error);
      curl_easy_setopt(easy, CURLOPT_USERAGENT, user_agent.c_str());
      curl_easy_setopt(easy, CURLOPT_NOSIGNAL, 1L); // thread-safe DNS timeouts
      curl_easy_setopt(easy, CURLOPT_FAILONERROR, 1L); // HTTP 4xx/5xx are errors, not documents
      curl_easy_setopt(easy, CURLOPT_FOLLOWLOCATION, 1L);
      curl_easy_setopt(easy, CURLOPT_MAXREDIRS, 20L);
      curl_easy_setopt(easy, CURLOPT_CONNECTTIMEOUT, 60L);
      // a stalled transfer fails instead of hanging: less than 1 byte/s for two minutes
      curl_easy_setopt(easy, CURLOPT_LOW_SPEED_LIMIT, 1L);
      curl_easy_setopt(easy, CURLOPT_LOW_SPEED_TIME, 120L);
#if LIBCURL_VERSION_NUM >= 0x075500 // 7.85.0
      static const std::string protocols = supportedProtocols({"http", "https", "ftp", "file"});
      static const std::string redirect_protocols = supportedProtocols({"http", "https", "ftp"}); // never to a local file
      ok = ok && curl_easy_setopt(easy, CURLOPT_PROTOCOLS_STR, protocols.c_str()) == CURLE_OK
              && curl_easy_setopt(easy, CURLOPT_REDIR_PROTOCOLS_STR, redirect_protocols.c_str()) == CURLE_OK;
#else
      ok = ok && curl_easy_setopt(easy, CURLOPT_PROTOCOLS, static_cast<long>(CURLPROTO_HTTP | CURLPROTO_HTTPS | CURLPROTO_FTP | CURLPROTO_FILE)) == CURLE_OK
              && curl_easy_setopt(easy, CURLOPT_REDIR_PROTOCOLS, static_cast<long>(CURLPROTO_HTTP | CURLPROTO_HTTPS | CURLPROTO_FTP)) == CURLE_OK;
#endif
      return ok;
    }

    // Throws the exception for the failed transfer @p result of @p easy from @p url (messages show neither its credentials
    // nor its query)
    [[noreturn]] void throwTransferError(CURL* easy, CURLcode result, const char* error, const std::string& url)
    {
      long response = 0;
      curl_easy_getinfo(easy, CURLINFO_RESPONSE_CODE, &response);
      const std::string name = HTTPInputSource::redactURL(url);
      const std::string message = HTTPInputSource::redactURLs(error[0] != '\0' ? std::string(error) : std::string(curl_easy_strerror(result)), url);
      if (response == 404 || response == 410 || result == CURLE_REMOTE_FILE_NOT_FOUND || result == CURLE_FILE_COULDNT_READ_FILE)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name, message);
      }
      // the message is logged, as callers may rethrow the exception as another one (e.g. MzMLFile::load(): a ParseError)
      OPENMS_LOG_ERROR << "Error reading " << name << ": " << message << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }

    // The value of one header of the last response of a transfer (a redirect that libcurl follows is a response too)
    struct ResponseHeader
    {
      std::string_view name; ///< in lower case, with the colon, e.g. "content-range:"
      std::string value;     ///< as sent: after the colon, with the line break

      // Collects the header from the responses to @p easy
      void collectFrom(CURL* easy)
      {
        curl_easy_setopt(easy, CURLOPT_HEADERFUNCTION, &ResponseHeader::header);
        curl_easy_setopt(easy, CURLOPT_HEADERDATA, this);
      }

      // The value without whitespace, in lower case (e.g. of Content-Encoding)
      std::string token() const
      {
        std::string result;
        for (const char c : value)
        {
          if (!std::isspace(static_cast<unsigned char>(c))) { result += static_cast<char>(std::tolower(static_cast<unsigned char>(c))); }
        }
        return result;
      }

      static size_t header(char* data, size_t size, size_t count, void* user)
      {
        auto& self = *static_cast<ResponseHeader*>(user);
        const std::string_view line(data, size * count);
        if (line.substr(0, 5) == "HTTP/") { self.value.clear(); } // the status line of a response
        else if (line.size() > self.name.size()
                 && std::equal(self.name.begin(), self.name.end(), line.begin(), [](char a, char b) { return a == std::tolower(static_cast<unsigned char>(b)); }))
        {
          self.value = line.substr(self.name.size());
        }
        return size * count;
      }
    };

    // A transfer whose received bytes are handed out on demand: readBytes() drives libcurl (multi interface) until data
    // arrives, so at most what one call of curl_multi_perform() receives is buffered
    class CurlInputStream :
      public InputStream
    {
public:
      explicit CurlInputStream(const std::string& url) :
        url_(url),
        name_(HTTPInputSource::redactURL(url))
      {
        CurlInit::ensure();
        easy_ = curl_easy_init();
        multi_ = curl_multi_init();
        if (easy_ == nullptr || multi_ == nullptr)
        {
          cleanup_();
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name_);
        }
        bool ok = setUpTransfer(easy_, url_, error_);
        curl_easy_setopt(easy_, CURLOPT_WRITEFUNCTION, write_);
        curl_easy_setopt(easy_, CURLOPT_WRITEDATA, this);
        // gzip transfer compression arrives undecoded and is decompressed as the stream is read (makeStream()).
        // libcurl would decode everything that one curl_multi_perform() receives at once, which for highly compressed
        // content is far more than one transfer step.
        headers_ = curl_slist_append(nullptr, "Accept-Encoding: gzip");
        ok = ok && headers_ != nullptr && curl_easy_setopt(easy_, CURLOPT_HTTPHEADER, headers_) == CURLE_OK;
        curl_easy_setopt(easy_, CURLOPT_HTTP_CONTENT_DECODING, 0L);
        content_encoding_.collectFrom(easy_);
        if (!ok || curl_multi_add_handle(multi_, easy_) != CURLM_OK)
        {
          OPENMS_LOG_ERROR << "Cannot set up the transfer of " << name_ << std::endl;
          cleanup_();
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name_);
        }
        added_ = true;
      }

      ~CurlInputStream() override
      {
        cleanup_();
      }

      UInt64 curPos() const override
      {
        return position_;
      }

      Size readBytes(unsigned char* to_fill, Size max_to_read) override
      {
        if (max_to_read == 0) { return 0; }
        while (pending_() == 0 && !done_) { pump_(); }
        if (pending_() == 0)
        {
          checkResult_();
          return 0;
        }
        const Size count = std::min(max_to_read, pending_());
        std::memcpy(to_fill, buffer_.data() + offset_, count);
        offset_ += count;
        position_ += count;
        if (offset_ == buffer_.size())
        {
          buffer_.clear();
          offset_ = 0;
        }
        return count;
      }

      void cancel() override
      {
        cancelled_ = true;
#if LIBCURL_VERSION_NUM >= 0x074400 // 7.68.0: wake up a curl_multi_poll() on the reading thread
        curl_multi_wakeup(multi_);
#endif
      }

      /// The first (at most) @p count unread bytes, without consuming them; throws if the transfer failed
      std::string_view peek(Size count)
      {
        while (pending_() < count && !done_) { pump_(); }
        if (pending_() == 0) { checkResult_(); }
        return std::string_view(buffer_.data() + offset_, std::min(count, pending_()));
      }

      /// Whether the response has gzip transfer compression (Content-Encoding); known after peek()
      bool transferCompressed() const
      {
        const std::string coding = content_encoding_.token();
        return coding == "gzip" || coding == "x-gzip";
      }

      CurlInputStream(const CurlInputStream&) = delete;
      CurlInputStream& operator=(const CurlInputStream&) = delete;

private:
      static size_t write_(char* data, size_t size, size_t count, void* user)
      {
        auto& self = *static_cast<CurlInputStream*>(user);
        if (self.offset_ > 0 && self.offset_ >= self.buffer_.size() / 2)
        {
          self.buffer_.erase(0, self.offset_);
          self.offset_ = 0;
        }
        self.buffer_.append(data, size * count);
        self.received_ += size * count;
        return size * count;
      }

      Size pending_() const
      {
        return buffer_.size() - offset_;
      }

      // Lets libcurl transfer what it can without blocking; if nothing arrived, waits a while for the connection
      void pump_()
      {
        if (cancelled_)
        {
          done_ = true; // the stream is not read after cancel()
          return;
        }
        const UInt64 received = received_;
        int running = 0;
        CURLMcode code = curl_multi_perform(multi_, &running);
        if (code != CURLM_OK)
        {
          OPENMS_LOG_ERROR << "Error reading " << name_ << ": " << curl_multi_strerror(code) << std::endl;
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name_);
        }
        if (running == 0)
        {
          int left = 0;
          while (CURLMsg* message = curl_multi_info_read(multi_, &left))
          {
            if (message->msg == CURLMSG_DONE) { result_ = message->data.result; }
          }
          done_ = true;
          return;
        }
        if (received_ == received)
        {
          // without curl_multi_wakeup(), a cancel() is noticed after the timeout
          code = curl_multi_poll(multi_, nullptr, 0, LIBCURL_VERSION_NUM >= 0x074400 ? 1000 : 200, nullptr);
          if (code != CURLM_OK)
          {
            OPENMS_LOG_ERROR << "Error reading " << name_ << ": " << curl_multi_strerror(code) << std::endl;
            throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name_);
          }
        }
      }

      // Throws if the (finished) transfer failed
      void checkResult_() const
      {
        if (result_ != CURLE_OK) { throwTransferError(easy_, result_, error_, url_); }
      }

      void cleanup_()
      {
        if (added_) { curl_multi_remove_handle(multi_, easy_); }
        if (multi_ != nullptr) { curl_multi_cleanup(multi_); }
        if (easy_ != nullptr) { curl_easy_cleanup(easy_); }
        curl_slist_free_all(headers_);
        added_ = false;
        multi_ = nullptr;
        easy_ = nullptr;
        headers_ = nullptr;
      }

      std::string url_;
      std::string name_;        ///< url_ for messages, without credentials (HTTPInputSource::redactURL())
      CURL* easy_ = nullptr;
      CURLM* multi_ = nullptr;
      curl_slist* headers_ = nullptr;
      ResponseHeader content_encoding_{"content-encoding:", {}};
      bool added_ = false;
      char error_[CURL_ERROR_SIZE] = {};
      std::string buffer_;      ///< received, not yet read bytes from offset_ on
      Size offset_ = 0;
      bool done_ = false;       ///< the transfer has ended
      CURLcode result_ = CURLE_OK;
      UInt64 position_ = 0;
      UInt64 received_ = 0;     ///< bytes received in total
      std::atomic<bool> cancelled_{false};
    };

    // A stream whose first bytes can be looked at (peek()) before they are read, for their magic bytes
    class PeekingInputStream :
      public InputStream
    {
public:
      explicit PeekingInputStream(std::unique_ptr<InputStream> source) :
        source_(std::move(source))
      {
      }

      UInt64 curPos() const override
      {
        return position_;
      }

      Size readBytes(unsigned char* to_fill, Size max_to_read) override
      {
        Size count = 0;
        if (offset_ < head_.size())
        {
          count = std::min(max_to_read, head_.size() - offset_);
          std::memcpy(to_fill, head_.data() + offset_, count);
          offset_ += count;
        }
        else
        {
          count = source_->readBytes(to_fill, max_to_read);
        }
        position_ += count;
        return count;
      }

      void cancel() override
      {
        source_->cancel();
      }

      /// The first (at most) @p count bytes, without consuming them; before the first readBytes()
      std::string_view peek(Size count)
      {
        unsigned char buffer[16];
        while (head_.size() < count)
        {
          const Size read = source_->readBytes(buffer, std::min<Size>(sizeof(buffer), count - head_.size()));
          if (read == 0) { break; }
          head_.append(reinterpret_cast<const char*>(buffer), read);
        }
        return std::string_view(head_).substr(0, count);
      }

private:
      std::unique_ptr<InputStream> source_;
      std::string head_; ///< bytes of source_ read by peek()
      Size offset_ = 0;  ///< of the next byte of head_ to read
      UInt64 position_ = 0;
    };

    // The compression format whose magic bytes start @p head, if any
    std::optional<DecompressingInputStream::Format> compression(std::string_view head)
    {
      return DecompressingInputStream::detect(reinterpret_cast<const unsigned char*>(head.data()), head.size());
    }

    // @p stream (which has peek()), decompressed if its content is gzip or bzip2 compressed (by the magic bytes)
    template <typename Stream>
    std::unique_ptr<InputStream> decompressed(std::unique_ptr<Stream> stream, const std::string& name)
    {
      const auto format = compression(stream->peek(3));
      if (format) { return std::make_unique<DecompressingInputStream>(std::move(stream), *format, name); }
      return stream;
    }

    // Collects @p limit bytes of a response body, from @p offset of the resource
    struct RangeCollector
    {
      CURL* easy = nullptr;
      std::string bytes;
      Int64 offset = 0;
      Int64 limit = 0;
      Int64 skip = 0;
      bool started = false;

      static size_t write(char* data, size_t size, size_t count, void* user)
      {
        auto& self = *static_cast<RangeCollector*>(user);
        if (!self.started)
        {
          // an HTTP server that ignores the range answers 200 with the whole resource instead of 206: skip to the offset
          self.started = true;
          long response = 0;
          const char* url = nullptr;
          curl_easy_getinfo(self.easy, CURLINFO_RESPONSE_CODE, &response);
          curl_easy_getinfo(self.easy, CURLINFO_EFFECTIVE_URL, &url);
          if (response == 200 && url != nullptr && StringUtils::hasPrefix(std::string(url), "http")) { self.skip = self.offset; }
        }
        const Int64 total = static_cast<Int64>(size * count);
        const Int64 begin = std::min(self.skip, total);
        self.skip -= begin;
        const Int64 wanted = std::min(total - begin, self.limit - static_cast<Int64>(self.bytes.size()));
        self.bytes.append(data + begin, static_cast<size_t>(wanted));
        // all bytes received: stop the transfer (only a server that ignores the range sends more)
        return static_cast<Int64>(self.bytes.size()) < self.limit ? size * count : 0;
      }
    };

    // A single blocking request (no stream), for objectSize() and readRange()
    class Request
    {
public:
      explicit Request(const std::string& url) :
        url_(url),
        name_(HTTPInputSource::redactURL(url))
      {
        CurlInit::ensure();
        easy_ = curl_easy_init();
        if (easy_ == nullptr || !setUpTransfer(easy_, url, error_))
        {
          OPENMS_LOG_ERROR << "Cannot set up the transfer of " << name_ << std::endl;
          if (easy_ != nullptr) { curl_easy_cleanup(easy_); }
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name_);
        }
      }

      ~Request()
      {
        curl_easy_cleanup(easy_);
      }

      Request(const Request&) = delete;
      Request& operator=(const Request&) = delete;

      CURL* easy() const
      {
        return easy_;
      }

      // Performs the request; @p accept tells whether a failed result is fine nevertheless
      template <typename Accept>
      void perform(Accept accept)
      {
        const CURLcode result = curl_easy_perform(easy_);
        if (result != CURLE_OK && !accept(result)) { throwTransferError(easy_, result, error_, url_); }
      }

      const std::string& url() const
      {
        return url_;
      }

private:
      std::string url_;
      std::string name_;
      CURL* easy_ = nullptr;
      char error_[CURL_ERROR_SIZE] = {};
    };

    // Receives the response to a request of the first byte of a resource (objectSize()): its Content-Range, not its body
    struct SizeProbe
    {
      ResponseHeader content_range{"content-range:", {}};

      // The body is not needed: the transfer stops at its first byte (CURLE_WRITE_ERROR). libcurl does not hand out the
      // body of a redirect it follows.
      static size_t write(char*, size_t, size_t, void*)
      {
        return 0;
      }

      // The size of the resource in the Content-Range ("bytes 0-0/<size>"; "bytes */<size>" for a range beyond its end);
      // -1 if there is none, or it is unknown ("*")
      Int64 total() const
      {
        const std::string& range = content_range.value;
        const auto slash = range.rfind('/');
        if (slash == std::string::npos) { return -1; }
        const char* begin = range.c_str() + slash + 1;
        char* end = nullptr;
        const long long size = std::strtoll(begin, &end, 10);
        return end != begin && size >= 0 ? static_cast<Int64>(size) : -1;
      }
    };

    // The size of the resource at @p url as reported for a request without body (FTP: SIZE; file: the file size); -1 if none
    curl_off_t sizeWithoutBody(const std::string& url)
    {
      Request request(url);
      curl_easy_setopt(request.easy(), CURLOPT_NOBODY, 1L);
      request.perform([](CURLcode) { return false; });
      curl_off_t size = -1;
      curl_easy_getinfo(request.easy(), CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &size);
      return size;
    }

    // The size of the HTTP(S) resource at @p url from the response to a request of its first byte; -1 if none. A GET, not
    // a HEAD request: a presigned URL (e.g. S3, GCS) signs the request method, so a URL presigned for GET is refused with HEAD.
    curl_off_t sizeFromFirstByte(const std::string& url)
    {
      Request request(url);
      SizeProbe probe;
      curl_easy_setopt(request.easy(), CURLOPT_RANGE, "0-0");
      curl_easy_setopt(request.easy(), CURLOPT_FAILONERROR, 0L); // the status is checked below: a 416 has the size too
      probe.content_range.collectFrom(request.easy());
      curl_easy_setopt(request.easy(), CURLOPT_WRITEFUNCTION, &SizeProbe::write);
      request.perform([](CURLcode result) { return result == CURLE_WRITE_ERROR; }); // stopped at the first byte
      const char* effective_url = nullptr;
      curl_easy_getinfo(request.easy(), CURLINFO_EFFECTIVE_URL, &effective_url);
      if (effective_url == nullptr || !StringUtils::hasPrefix(std::string(effective_url), "http"))
      {
        return sizeWithoutBody(url); // redirected to FTP
      }
      long response = 0;
      curl_easy_getinfo(request.easy(), CURLINFO_RESPONSE_CODE, &response);
      if (response >= 400 && response != 416)
      {
        const std::string error = "The requested URL returned error: " + std::to_string(response);
        throwTransferError(request.easy(), CURLE_HTTP_RETURNED_ERROR, error.c_str(), request.url());
      }
      curl_off_t size = -1;
      if (response == 206 || response == 416) // 416: the resource is empty
      {
        size = probe.total();
      }
      else if (response == 200) // a server that ignores the range sends the whole resource
      {
        curl_easy_getinfo(request.easy(), CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &size);
      }
      return size;
    }
  } // namespace

  HTTPInputSource::HTTPInputSource(const std::string& url) :
    url_(url)
  {
    if (!isURL(url) && !StringUtils::hasPrefix(url, "file://"))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not an HTTP, HTTPS, FTP or file URL: '" + redactURL(url) + "'");
    }
  }

  std::unique_ptr<InputStream> HTTPInputSource::makeStream() const
  {
    const std::string name = redactURL(url_);
    auto stream = std::make_unique<CurlInputStream>(url_);
    const auto format = compression(stream->peek(3)); // the response headers are in with its first bytes
    if (stream->transferCompressed() && format == DecompressingInputStream::Format::GZIP)
    {
      // gzip transfer compression, around the resource, which may be compressed itself (e.g. .mzML.gz)
      auto transfer = std::make_unique<DecompressingInputStream>(std::move(stream), DecompressingInputStream::Format::GZIP, name);
      return decompressed(std::make_unique<PeekingInputStream>(std::move(transfer)), name);
    }
    // a gzip or bzip2 compressed resource (e.g. .mzML.gz), as for a local file
    return decompressed(std::move(stream), name);
  }

  Int64 HTTPInputSource::objectSize(const std::string& url)
  {
    const HTTPInputSource source(url);
    const curl_off_t size = StringUtils::hasPrefix(url, "http") ? sizeFromFirstByte(url) : sizeWithoutBody(url);
    if (size < 0)
    {
      const std::string name = redactURL(url);
      OPENMS_LOG_ERROR << "Error reading " << name << ": the server does not report its size" << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
    }
    return static_cast<Int64>(size);
  }

  std::string HTTPInputSource::readRange(const std::string& url, Int64 offset, Int64 length)
  {
    const HTTPInputSource source(url);
    if (length <= 0) { return {}; }
    Request request(url);
    RangeCollector collector;
    collector.easy = request.easy();
    collector.offset = offset;
    collector.limit = length;
    const std::string range = std::to_string(offset) + "-" + std::to_string(offset + length - 1);
    curl_easy_setopt(request.easy(), CURLOPT_RANGE, range.c_str());
    curl_easy_setopt(request.easy(), CURLOPT_WRITEFUNCTION, &RangeCollector::write);
    curl_easy_setopt(request.easy(), CURLOPT_WRITEDATA, &collector);
    // CURLE_WRITE_ERROR: RangeCollector stopped the transfer once it had all bytes
    request.perform([&](CURLcode result) { return result == CURLE_WRITE_ERROR && static_cast<Int64>(collector.bytes.size()) == length; });
    return collector.bytes;
  }

  bool HTTPInputSource::isURL(const std::string& url)
  {
    return StringUtils::hasPrefix(url, "http://") || StringUtils::hasPrefix(url, "https://") || StringUtils::hasPrefix(url, "ftp://");
  }

  std::string HTTPInputSource::redactURL(const std::string& url)
  {
    // a URL starts with a scheme (RFC 3986: a letter, then letters, digits, '+', '-' or '.') and "://"; at least two
    // characters, so that a Windows drive ("C://data/...") stays a file name
    const auto scheme_end = url.find("://");
    const auto scheme_char = [](unsigned char c) { return std::isalnum(c) || c == '+' || c == '-' || c == '.'; };
    if (scheme_end == std::string::npos || scheme_end < 2 || !std::isalpha(static_cast<unsigned char>(url[0]))
        || !std::all_of(url.begin(), url.begin() + scheme_end, scheme_char))
    {
      return url;
    }
    const Size authority_begin = scheme_end + 3;
    const Size authority_end = std::min(url.find_first_of("/?#", authority_begin), url.size());
    const std::string authority = url.substr(authority_begin, authority_end - authority_begin);
    const auto at = authority.rfind('@');
    const Size path_end = std::min(url.find_first_of("?#", authority_end), url.size());
    return url.substr(0, authority_begin) + (at == std::string::npos ? authority : authority.substr(at + 1))
           + url.substr(authority_end, path_end - authority_end);
  }

  std::string HTTPInputSource::redactURLs(const std::string& text)
  {
    const auto scheme_char = [](unsigned char c) { return std::isalnum(c) || c == '+' || c == '-' || c == '.'; };
    const auto url_end = [](unsigned char c) { return std::isspace(c) || c == '\'' || c == '"' || c == '`' || c == '<' || c == '>'; };
    std::string result;
    Size copied = 0; // text before this is in result
    for (Size separator = text.find("://"); separator != std::string::npos; separator = text.find("://", separator + 3))
    {
      Size begin = separator;
      while (begin > copied && scheme_char(static_cast<unsigned char>(text[begin - 1]))) { --begin; }
      if (begin == separator) { continue; } // no scheme
      Size end = separator + 3;
      while (end < text.size() && !url_end(static_cast<unsigned char>(text[end]))) { ++end; }
      result += text.substr(copied, begin - copied) + redactURL(text.substr(begin, end - begin));
      copied = end;
      separator = end - 3; // the next search starts after this URL
    }
    return result + text.substr(copied);
  }

  std::string HTTPInputSource::redactURLs(const std::string& text, const std::string& url)
  {
    std::string result = redactURLs(text);
    if (redactURL(url) == url) { return result; } // a local file name, or nothing to remove
    // the parts that redactURL() removes: the user name and password, each query value, the fragment; as given and decoded
    const Size authority_begin = url.find("://") + 3;
    const Size authority_end = std::min(url.find_first_of("/?#", authority_begin), url.size());
    const std::string_view authority = std::string_view(url).substr(authority_begin, authority_end - authority_begin);
    std::vector<std::string_view> parts;
    const auto at = authority.rfind('@');
    if (at != std::string_view::npos)
    {
      const std::string_view user_info = authority.substr(0, at);
      const auto colon = user_info.find(':');
      parts.push_back(user_info.substr(0, colon));
      if (colon != std::string_view::npos) { parts.push_back(user_info.substr(colon + 1)); }
    }
    const Size query_begin = std::min(url.find_first_of("?#", authority_end), url.size());
    const Size fragment_begin = std::min(url.find('#', query_begin), url.size());
    std::string_view query = std::string_view(url).substr(query_begin, fragment_begin - query_begin);
    while (!query.empty())
    {
      query.remove_prefix(1); // '?' or '&'
      const std::string_view item = query.substr(0, query.find('&'));
      const auto equals = item.find('=');
      if (equals != std::string_view::npos) { parts.push_back(item.substr(equals + 1)); }
      query.remove_prefix(item.size());
    }
    if (fragment_begin + 1 < url.size()) { parts.push_back(std::string_view(url).substr(fragment_begin + 1)); }
    std::vector<std::string> secrets;
    for (const std::string_view part : parts)
    {
      std::string plus_as_space(part);
      std::replace(plus_as_space.begin(), plus_as_space.end(), '+', ' ');
      for (std::string secret : {std::string(part), percentDecoded(part), percentDecoded(plus_as_space)})
      {
        if (!secret.empty() && std::find(secrets.begin(), secrets.end(), secret) == secrets.end()) { secrets.push_back(std::move(secret)); }
      }
    }
    // longer ones first, so that a part of one does not break up another
    std::sort(secrets.begin(), secrets.end(), [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
    const auto word_char = [](unsigned char c) { return std::isalnum(c) != 0; };
    for (const std::string& secret : secrets)
    {
      Size found = result.find(secret);
      while (found != std::string::npos)
      {
        const Size end = found + secret.size();
        // only a whole word, not a part of a longer one (e.g. a value "1" in "port 1234")
        const bool part_of_word = (found > 0 && word_char(result[found - 1]) && word_char(secret.front()))
                                  || (end < result.size() && word_char(result[end]) && word_char(secret.back()));
        if (part_of_word)
        {
          found = result.find(secret, found + 1);
          continue;
        }
        result.replace(found, secret.size(), "***");
        found = result.find(secret, found + 3);
      }
    }
    return result;
  }

} // namespace OpenMS
