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
#include <cstring>
#include <initializer_list>
#include <string_view>

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

    // A transfer whose received bytes are handed out on demand: readBytes() drives libcurl (multi interface) until data
    // arrives, so at most what one call of curl_multi_perform() receives is buffered
    class CurlInputStream :
      public InputStream
    {
public:
      explicit CurlInputStream(const std::string& url) :
        url_(url)
      {
        CurlInit::ensure();
        easy_ = curl_easy_init();
        multi_ = curl_multi_init();
        if (easy_ == nullptr || multi_ == nullptr)
        {
          cleanup_();
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
        }
        user_agent_ = "OpenMS/" + VersionInfo::getVersion();
        // without the URL or the protocol restriction, the transfer must not start
        bool ok = curl_easy_setopt(easy_, CURLOPT_URL, url_.c_str()) == CURLE_OK;
        curl_easy_setopt(easy_, CURLOPT_WRITEFUNCTION, write_);
        curl_easy_setopt(easy_, CURLOPT_WRITEDATA, this);
        curl_easy_setopt(easy_, CURLOPT_ERRORBUFFER, error_);
        curl_easy_setopt(easy_, CURLOPT_USERAGENT, user_agent_.c_str());
        curl_easy_setopt(easy_, CURLOPT_NOSIGNAL, 1L); // thread-safe DNS timeouts
        curl_easy_setopt(easy_, CURLOPT_FAILONERROR, 1L); // HTTP 4xx/5xx are errors, not documents
        curl_easy_setopt(easy_, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(easy_, CURLOPT_MAXREDIRS, 20L);
        curl_easy_setopt(easy_, CURLOPT_ACCEPT_ENCODING, ""); // any transfer compression libcurl can decode
        curl_easy_setopt(easy_, CURLOPT_CONNECTTIMEOUT, 60L);
        // a stalled transfer fails instead of hanging: less than 1 byte/s for two minutes
        curl_easy_setopt(easy_, CURLOPT_LOW_SPEED_LIMIT, 1L);
        curl_easy_setopt(easy_, CURLOPT_LOW_SPEED_TIME, 120L);
#if LIBCURL_VERSION_NUM >= 0x075500 // 7.85.0
        const std::string protocols = supportedProtocols({"http", "https", "ftp", "file"});
        const std::string redirect_protocols = supportedProtocols({"http", "https", "ftp"}); // never to a local file
        ok = ok && curl_easy_setopt(easy_, CURLOPT_PROTOCOLS_STR, protocols.c_str()) == CURLE_OK
                && curl_easy_setopt(easy_, CURLOPT_REDIR_PROTOCOLS_STR, redirect_protocols.c_str()) == CURLE_OK;
#else
        ok = ok && curl_easy_setopt(easy_, CURLOPT_PROTOCOLS, static_cast<long>(CURLPROTO_HTTP | CURLPROTO_HTTPS | CURLPROTO_FTP | CURLPROTO_FILE)) == CURLE_OK
                && curl_easy_setopt(easy_, CURLOPT_REDIR_PROTOCOLS, static_cast<long>(CURLPROTO_HTTP | CURLPROTO_HTTPS | CURLPROTO_FTP)) == CURLE_OK;
#endif
        if (!ok || curl_multi_add_handle(multi_, easy_) != CURLM_OK)
        {
          OPENMS_LOG_ERROR << "Cannot set up the transfer of " << url_ << std::endl;
          cleanup_();
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
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
          OPENMS_LOG_ERROR << "Error reading " << url_ << ": " << curl_multi_strerror(code) << std::endl;
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
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
            OPENMS_LOG_ERROR << "Error reading " << url_ << ": " << curl_multi_strerror(code) << std::endl;
            throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
          }
        }
      }

      // Throws if the (finished) transfer failed
      void checkResult_() const
      {
        if (result_ == CURLE_OK) { return; }
        long response = 0;
        curl_easy_getinfo(easy_, CURLINFO_RESPONSE_CODE, &response);
        const std::string message = error_[0] != '\0' ? std::string(error_) : std::string(curl_easy_strerror(result_));
        if (response == 404 || response == 410 || result_ == CURLE_REMOTE_FILE_NOT_FOUND || result_ == CURLE_FILE_COULDNT_READ_FILE)
        {
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_, message);
        }
        // the message is logged, as callers may rethrow the exception as another one (e.g. MzMLFile::load(): a ParseError)
        OPENMS_LOG_ERROR << "Error reading " << url_ << ": " << message << std::endl;
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, url_);
      }

      void cleanup_()
      {
        if (added_) { curl_multi_remove_handle(multi_, easy_); }
        if (multi_ != nullptr) { curl_multi_cleanup(multi_); }
        if (easy_ != nullptr) { curl_easy_cleanup(easy_); }
        added_ = false;
        multi_ = nullptr;
        easy_ = nullptr;
      }

      std::string url_;
      std::string user_agent_;
      CURL* easy_ = nullptr;
      CURLM* multi_ = nullptr;
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
  } // namespace

  HTTPInputSource::HTTPInputSource(const std::string& url) :
    url_(url)
  {
    if (!isURL(url) && !StringUtils::hasPrefix(url, "file://"))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not an HTTP, HTTPS, FTP or file URL: '" + url + "'");
    }
  }

  std::unique_ptr<InputStream> HTTPInputSource::makeStream() const
  {
    auto stream = std::make_unique<CurlInputStream>(url_);
    // gzip or bzip2 compressed: by the magic bytes
    const std::string_view head = stream->peek(3);
    const auto format = DecompressingInputStream::detect(reinterpret_cast<const unsigned char*>(head.data()), head.size());
    if (format)
    {
      return std::make_unique<DecompressingInputStream>(std::move(stream), *format, url_);
    }
    return stream;
  }

  bool HTTPInputSource::isURL(const std::string& url)
  {
    return StringUtils::hasPrefix(url, "http://") || StringUtils::hasPrefix(url, "https://") || StringUtils::hasPrefix(url, "ftp://");
  }

} // namespace OpenMS
