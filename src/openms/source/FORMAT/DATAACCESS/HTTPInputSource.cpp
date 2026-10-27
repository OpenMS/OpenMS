// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/DecompressingInputStream.h>
#include <OpenMS/SYSTEM/CurlInit.h>

#include <curl/curl.h>

#include <algorithm>
#include <cstring>

namespace OpenMS
{
  namespace
  {
    bool hasSchemePrefix(const std::string& url, const std::string& prefix)
    {
      return url.size() >= prefix.size() && StringUtils::toLowered(url.substr(0, prefix.size())) == prefix;
    }

    /// An easy handle with the options all requests share
    class Request
    {
public:
      explicit Request(const std::string& url) :
        url_(url)
      {
        CurlInit::ensure();
        handle_ = curl_easy_init();
        if (handle_ == nullptr)
        {
          OPENMS_LOG_ERROR << "Cannot read " << url_ << ": failed to initialize libcurl" << std::endl;
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
        }
        curl_easy_setopt(handle_, CURLOPT_URL, url_.c_str());
        curl_easy_setopt(handle_, CURLOPT_ERRORBUFFER, error_);
        curl_easy_setopt(handle_, CURLOPT_FOLLOWLOCATION, 1L);
        // HTTP status >= 400 is an error, not a body to parse
        curl_easy_setopt(handle_, CURLOPT_FAILONERROR, 1L);
        curl_easy_setopt(handle_, CURLOPT_NOSIGNAL, 1L);
        curl_easy_setopt(handle_, CURLOPT_USERAGENT, "OpenMS");
        curl_easy_setopt(handle_, CURLOPT_CONNECTTIMEOUT, 60L);
        // give up on a transfer that delivers nothing for five minutes
        curl_easy_setopt(handle_, CURLOPT_LOW_SPEED_LIMIT, 1L);
        curl_easy_setopt(handle_, CURLOPT_LOW_SPEED_TIME, 300L);
        // these protocols only, also after redirects (no file:// or others)
#if LIBCURL_VERSION_NUM >= 0x075500
        curl_easy_setopt(handle_, CURLOPT_PROTOCOLS_STR, "http,https,ftp");
        curl_easy_setopt(handle_, CURLOPT_REDIR_PROTOCOLS_STR, "http,https,ftp");
#else
        curl_easy_setopt(handle_, CURLOPT_PROTOCOLS, CURLPROTO_HTTP | CURLPROTO_HTTPS | CURLPROTO_FTP);
        curl_easy_setopt(handle_, CURLOPT_REDIR_PROTOCOLS, CURLPROTO_HTTP | CURLPROTO_HTTPS | CURLPROTO_FTP);
#endif
      }

      ~Request()
      {
        curl_easy_cleanup(handle_);
      }

      Request(const Request&) = delete;
      Request& operator=(const Request&) = delete;

      CURL* handle() const
      {
        return handle_;
      }

      /// Throws the exception for the failed transfer @p result
      [[noreturn]] void fail(CURLcode result) const
      {
        long status = 0;
        curl_easy_getinfo(handle_, CURLINFO_RESPONSE_CODE, &status);
        const std::string reason = error_[0] != '\0' ? std::string(error_) : std::string(curl_easy_strerror(result));
        // logged, as callers may rethrow it as another exception (e.g. MzMLFile::load(): a ParseError)
        OPENMS_LOG_ERROR << "Cannot read " << url_ << ": " << reason << std::endl;
        if (result == CURLE_REMOTE_FILE_NOT_FOUND || (result == CURLE_HTTP_RETURNED_ERROR && (status == 404 || status == 410)))
        {
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_, reason);
        }
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
      }

private:
      std::string url_;
      CURL* handle_ = nullptr;
      char error_[CURL_ERROR_SIZE] = {};
    };

    /**
      @brief The response body of a GET request, pulled by the reader

      curl's multi interface moves the transfer forward on the calling thread whenever the received bytes are used up,
      so the download proceeds as fast as the parser (or the read-ahead thread) reads.
    */
    class CurlInputStream :
      public InputStream
    {
public:
      explicit CurlInputStream(const std::string& url) :
        request_(url)
      {
        CURL* handle = request_.handle();
        curl_easy_setopt(handle, CURLOPT_WRITEFUNCTION, &CurlInputStream::write_);
        curl_easy_setopt(handle, CURLOPT_WRITEDATA, this);
        // a server may compress the transfer (Content-Encoding); the body arrives decoded
        curl_easy_setopt(handle, CURLOPT_ACCEPT_ENCODING, "");
        multi_ = curl_multi_init();
        if (multi_ == nullptr || curl_multi_add_handle(multi_, handle) != CURLM_OK)
        {
          if (multi_ != nullptr) { curl_multi_cleanup(multi_); }
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url);
        }
      }

      ~CurlInputStream() override
      {
        curl_multi_remove_handle(multi_, request_.handle());
        curl_multi_cleanup(multi_);
      }

      UInt64 curPos() const override
      {
        return position_;
      }

      Size readBytes(unsigned char* to_fill, Size max_to_read) override
      {
        if (!fill_(1)) { return 0; }
        const Size count = std::min(max_to_read, pending_.size() - offset_);
        std::memcpy(to_fill, pending_.data() + offset_, count);
        offset_ += count;
        position_ += count;
        return count;
      }

      /// Up to @p count of the next bytes, without consuming them; fewer only at the end of the body
      std::string peek(Size count)
      {
        fill_(count);
        return pending_.substr(offset_, count);
      }

private:
      static size_t write_(char* data, size_t size, size_t count, void* self)
      {
        static_cast<CurlInputStream*>(self)->pending_.append(data, size * count);
        return size * count;
      }

      /// Receives until @p count bytes are pending or the transfer ended; whether any byte is pending
      bool fill_(Size count)
      {
        if (offset_ == pending_.size())
        {
          pending_.clear();
          offset_ = 0;
        }
        while (pending_.size() - offset_ < count && !done_)
        {
          int running = 0;
          const CURLMcode status = curl_multi_perform(multi_, &running);
          if (status != CURLM_OK)
          {
            OPENMS_LOG_ERROR << "libcurl error: " << curl_multi_strerror(status) << std::endl;
            throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "");
          }
          if (running == 0)
          {
            done_ = true;
            int left = 0;
            while (const CURLMsg* message = curl_multi_info_read(multi_, &left))
            {
              if (message->msg == CURLMSG_DONE && message->data.result != CURLE_OK) { request_.fail(message->data.result); }
            }
          }
          else if (pending_.size() - offset_ < count)
          {
            curl_multi_poll(multi_, nullptr, 0, 1000, nullptr);
          }
        }
        return offset_ < pending_.size();
      }

      Request request_;
      CURLM* multi_ = nullptr;
      std::string pending_; ///< received bytes; those before offset_ are consumed
      Size offset_ = 0;
      bool done_ = false;
      UInt64 position_ = 0;
    };

    /// Collects @p limit bytes of a response body, from @p offset of the file
    struct RangeCollector
    {
      CURL* handle = nullptr;
      std::string bytes;
      Int64 offset = 0;
      Int64 limit = 0;
      Int64 skip = 0;
      bool started = false;

      static size_t write(char* data, size_t size, size_t count, void* self)
      {
        auto& collector = *static_cast<RangeCollector*>(self);
        if (!collector.started)
        {
          // an HTTP server that ignores the range answers 200 with the whole file instead of 206: skip to the offset
          collector.started = true;
          long status = 0;
          const char* url = nullptr;
          curl_easy_getinfo(collector.handle, CURLINFO_RESPONSE_CODE, &status);
          curl_easy_getinfo(collector.handle, CURLINFO_EFFECTIVE_URL, &url);
          if (status == 200 && url != nullptr && hasSchemePrefix(url, "http")) { collector.skip = collector.offset; }
        }
        const Int64 total = static_cast<Int64>(size * count);
        const Int64 begin = std::min(collector.skip, total);
        collector.skip -= begin;
        const Int64 wanted = std::min(total - begin, collector.limit - static_cast<Int64>(collector.bytes.size()));
        collector.bytes.append(data + begin, static_cast<size_t>(wanted));
        // all bytes received: stop the transfer (only a server ignoring the range sends more)
        return static_cast<Int64>(collector.bytes.size()) < collector.limit ? size * count : 0;
      }
    };
  } // namespace

  HTTPInputSource::HTTPInputSource(const std::string& url) :
    url_(url)
  {
    if (!isHTTPURI(url))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not an HTTP, HTTPS or FTP URL: '" + url + "'");
    }
  }

  HTTPInputSource::~HTTPInputSource() = default;

  std::unique_ptr<InputStream> HTTPInputSource::makeStream() const
  {
    auto stream = std::make_unique<CurlInputStream>(url_);
    // starts the download: a missing file throws here, not in the parser
    const std::string magic = stream->peek(3);
    const auto format = DecompressingInputStream::detect(reinterpret_cast<const unsigned char*>(magic.data()), magic.size());
    if (format)
    {
      return std::make_unique<DecompressingInputStream>(std::move(stream), *format, url_);
    }
    return stream;
  }

  bool HTTPInputSource::isHTTPURI(const std::string& url)
  {
    return hasSchemePrefix(url, "http://") || hasSchemePrefix(url, "https://") || hasSchemePrefix(url, "ftp://");
  }

  Int64 HTTPInputSource::objectSize(const std::string& url)
  {
    const HTTPInputSource source(url);
    Request request(url);
    curl_easy_setopt(request.handle(), CURLOPT_NOBODY, 1L);
    const CURLcode result = curl_easy_perform(request.handle());
    if (result != CURLE_OK) { request.fail(result); }
    curl_off_t size = -1;
    curl_easy_getinfo(request.handle(), CURLINFO_CONTENT_LENGTH_DOWNLOAD_T, &size);
    if (size < 0)
    {
      OPENMS_LOG_ERROR << "Cannot read " << url << ": the server does not report its size" << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url);
    }
    return static_cast<Int64>(size);
  }

  std::string HTTPInputSource::readRange(const std::string& url, Int64 offset, Int64 length)
  {
    const HTTPInputSource source(url);
    if (length <= 0) { return {}; }
    Request request(url);
    RangeCollector collector;
    collector.handle = request.handle();
    collector.offset = offset;
    collector.limit = length;
    const std::string range = std::to_string(offset) + "-" + std::to_string(offset + length - 1);
    curl_easy_setopt(request.handle(), CURLOPT_RANGE, range.c_str());
    curl_easy_setopt(request.handle(), CURLOPT_WRITEFUNCTION, &RangeCollector::write);
    curl_easy_setopt(request.handle(), CURLOPT_WRITEDATA, &collector);
    const CURLcode result = curl_easy_perform(request.handle());
    // CURLE_WRITE_ERROR: stopped by RangeCollector once all bytes were received
    if (result != CURLE_OK && !(result == CURLE_WRITE_ERROR && static_cast<Int64>(collector.bytes.size()) == length))
    {
      request.fail(result);
    }
    return collector.bytes;
  }

} // namespace OpenMS
