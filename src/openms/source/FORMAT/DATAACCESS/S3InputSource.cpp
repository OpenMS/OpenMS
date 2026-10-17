// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/S3InputSource.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/HANDLERS/StringManager.h>

#include <xercesc/util/BinInputStream.hpp>

#ifdef OPENMS_WITH_ARROW_S3
#include <arrow/buffer.h>
#include <arrow/filesystem/s3fs.h>
#include <arrow/io/buffered.h>
#include <arrow/io/interfaces.h>
#include <arrow/memory_pool.h>
#endif

#include <bzlib.h>
#include <zlib.h>

#include <cstdlib>
#include <mutex>
#include <string_view>

namespace OpenMS
{

  struct S3InputSource::Impl
  {
    std::string uri;
#ifdef OPENMS_WITH_ARROW_S3
    std::shared_ptr<arrow::io::RandomAccessFile> object;
#endif
  };

  namespace
  {
#ifdef OPENMS_WITH_ARROW_S3
    // Each read of an S3 object is a request: the streams read the object in blocks of this size
    constexpr int64_t READ_BLOCK_SIZE = 8 * 1024 * 1024;

    // Arrow's S3 APIs are initialized before their first use, and must be finalized before the process exits (it
    // crashes at shutdown otherwise, see finalizeS3AtExit())
    void ensureS3Initialized()
    {
      static std::once_flag once;
      static arrow::Status status;
      std::call_once(once, [] { status = arrow::fs::EnsureS3Initialized(); });
      if (!status.ok())
      {
        OPENMS_LOG_ERROR << "Cannot initialize S3: " << status.ToString() << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot initialize S3: " + status.ToString());
      }
    }

    // Finalizes S3 at exit. Called once an S3 filesystem (with its client) exists: the finalization uses statics that the
    // first client creates, and exit handlers and static destructors run in the reverse order of their registration,
    // so the handler is registered after them to run before they are destroyed.
    void finalizeS3AtExit()
    {
      static std::once_flag once;
      std::call_once(once, [] { std::atexit([] { (void)arrow::fs::EnsureS3Finalized(); }); });
    }

    // Opens the S3 object @p uri
    std::shared_ptr<arrow::io::RandomAccessFile> openObject(const std::string& uri)
    {
      ensureS3Initialized();
      std::string path;
      auto options = arrow::fs::S3Options::FromUri(uri, &path);
      // the errors are logged, as callers may rethrow them as other exceptions (e.g. MzMLFile::load(): a ParseError)
      if (!options.ok())
      {
        OPENMS_LOG_ERROR << "Invalid S3 URI " << uri << ": " << options.status().ToString() << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Invalid S3 URI '" + uri + "': " + options.status().ToString());
      }
      auto filesystem = arrow::fs::S3FileSystem::Make(*options);
      if (!filesystem.ok())
      {
        OPENMS_LOG_ERROR << "Cannot open " << uri << ": " << filesystem.status().ToString() << std::endl;
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, uri, filesystem.status().ToString());
      }
      finalizeS3AtExit();
      auto object = (*filesystem)->OpenInputFile(path);
      if (!object.ok())
      {
        OPENMS_LOG_ERROR << "Cannot open " << uri << ": " << object.status().ToString() << std::endl;
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, uri, object.status().ToString());
      }
      return *object;
    }

    // The bytes of an Arrow stream, for Xerces
    class ArrowBinInputStream :
      public xercesc::BinInputStream
    {
public:
      explicit ArrowBinInputStream(std::shared_ptr<arrow::io::InputStream> stream) :
        stream_(std::move(stream))
      {
      }

      XMLFilePos curPos() const override
      {
        return position_;
      }

      XMLSize_t readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read) override
      {
        auto bytes_read = stream_->Read(static_cast<int64_t>(max_to_read), to_fill);
        if (!bytes_read.ok())
        {
          OPENMS_LOG_ERROR << "Error reading from S3: " << bytes_read.status().ToString() << std::endl;
          return 0;
        }
        position_ += *bytes_read;
        return static_cast<XMLSize_t>(*bytes_read);
      }

      const XMLCh* getContentType() const override
      {
        return nullptr;
      }

private:
      std::shared_ptr<arrow::io::InputStream> stream_;
      XMLFilePos position_ = 0;
    };

    // The decompressed bytes of a gzip (@p GZIP true) or bzip2 compressed Arrow stream, for Xerces
    class DecompressingBinInputStream :
      public xercesc::BinInputStream
    {
public:
      DecompressingBinInputStream(std::shared_ptr<arrow::io::InputStream> stream, bool gzip) :
        stream_(std::move(stream)),
        gzip_(gzip)
      {
        // windowBits 15 + 16: gzip format
        if (gzip_ ? inflateInit2(&z_stream_, 15 + 16) != Z_OK : BZ2_bzDecompressInit(&bz_stream_, 0, 0) != BZ_OK)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           std::string("Cannot initialize the ") + (gzip_ ? "gzip" : "bzip2") + " decompression");
        }
      }

      ~DecompressingBinInputStream() override
      {
        if (gzip_)
        {
          inflateEnd(&z_stream_);
        }
        else
        {
          BZ2_bzDecompressEnd(&bz_stream_);
        }
      }

      XMLFilePos curPos() const override
      {
        return position_;
      }

      XMLSize_t readBytes(XMLByte* const to_fill, const XMLSize_t max_to_read) override
      {
        XMLSize_t written = 0;
        while (written < max_to_read && !finished_)
        {
          if (available_ == 0)
          {
            auto bytes_read = stream_->Read(sizeof(buffer_), buffer_);
            if (!bytes_read.ok())
            {
              OPENMS_LOG_ERROR << "Error reading from S3: " << bytes_read.status().ToString() << std::endl;
              finished_ = true;
              break;
            }
            if (*bytes_read == 0) // truncated stream
            {
              finished_ = true;
              break;
            }
            next_ = buffer_;
            available_ = static_cast<size_t>(*bytes_read);
          }
          const size_t space = max_to_read - written;
          size_t consumed = 0;
          size_t produced = 0;
          bool ok = true;
          if (gzip_)
          {
            z_stream_.next_in = reinterpret_cast<Bytef*>(next_);
            z_stream_.avail_in = static_cast<uInt>(available_);
            z_stream_.next_out = reinterpret_cast<Bytef*>(to_fill + written);
            z_stream_.avail_out = static_cast<uInt>(space);
            const int ret = inflate(&z_stream_, Z_NO_FLUSH);
            consumed = available_ - z_stream_.avail_in;
            produced = space - z_stream_.avail_out;
            finished_ = ret == Z_STREAM_END;
            ok = ret == Z_OK || ret == Z_STREAM_END || ret == Z_BUF_ERROR;
          }
          else
          {
            bz_stream_.next_in = next_;
            bz_stream_.avail_in = static_cast<unsigned int>(available_);
            bz_stream_.next_out = reinterpret_cast<char*>(to_fill + written);
            bz_stream_.avail_out = static_cast<unsigned int>(space);
            const int ret = BZ2_bzDecompress(&bz_stream_);
            consumed = available_ - bz_stream_.avail_in;
            produced = space - bz_stream_.avail_out;
            finished_ = ret == BZ_STREAM_END;
            ok = ret == BZ_OK || ret == BZ_STREAM_END;
          }
          if (!ok || (consumed == 0 && produced == 0 && !finished_)) // an error, or no progress
          {
            OPENMS_LOG_ERROR << "Error decompressing the " << (gzip_ ? "gzip" : "bzip2") << " data read from S3" << std::endl;
            finished_ = true;
          }
          next_ += consumed;
          available_ -= consumed;
          written += produced;
        }
        position_ += written;
        return written;
      }

      const XMLCh* getContentType() const override
      {
        return nullptr;
      }

private:
      std::shared_ptr<arrow::io::InputStream> stream_;
      bool gzip_;
      z_stream z_stream_{};
      bz_stream bz_stream_{};
      char buffer_[64 * 1024];
      char* next_ = buffer_;   ///< the next compressed byte in buffer_
      size_t available_ = 0;   ///< the number of compressed bytes from next_ on
      bool finished_ = false;
      XMLFilePos position_ = 0;
    };
#endif
  } // namespace

  S3InputSource::S3InputSource(const std::string& uri) :
    impl_(std::make_unique<Impl>())
  {
    if (!isS3URI(uri))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not an S3 URI (s3://bucket/key): '" + uri + "'");
    }
    if (!isAvailable())
    {
      const std::string message = "Cannot read '" + uri + "': OpenMS was built without S3 support (its Arrow library has none, "
                                  "e.g. vcpkg's arrow without the feature 's3').";
      OPENMS_LOG_ERROR << message << std::endl; // callers may rethrow it as another exception (e.g. MzMLFile::load())
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }
    impl_->uri = uri;
#ifdef OPENMS_WITH_ARROW_S3
    impl_->object = openObject(uri);
#endif
    setSystemId(Internal::StringManager::convert(uri).c_str());
  }

  S3InputSource::~S3InputSource() = default;

  xercesc::BinInputStream* S3InputSource::makeStream() const
  {
#ifdef OPENMS_WITH_ARROW_S3
    // a stream of its own for each call, read in large blocks
    auto size = impl_->object->GetSize();
    if (!size.ok())
    {
      OPENMS_LOG_ERROR << "Error reading " << impl_->uri << ": " << size.status().ToString() << std::endl;
      return nullptr;
    }
    auto raw = arrow::io::RandomAccessFile::GetStream(impl_->object, 0, *size);
    if (!raw.ok())
    {
      OPENMS_LOG_ERROR << "Error reading " << impl_->uri << ": " << raw.status().ToString() << std::endl;
      return nullptr;
    }
    auto stream = arrow::io::BufferedInputStream::Create(READ_BLOCK_SIZE, arrow::default_memory_pool(), *raw);
    if (!stream.ok())
    {
      OPENMS_LOG_ERROR << "Error reading " << impl_->uri << ": " << stream.status().ToString() << std::endl;
      return nullptr;
    }
    // gzip or bzip2 compressed: by the magic bytes
    auto head = (*stream)->Peek(3);
    const std::string_view magic = head.ok() ? *head : std::string_view();
    if (magic.size() >= 2 && static_cast<unsigned char>(magic[0]) == 0x1F && static_cast<unsigned char>(magic[1]) == 0x8B)
    {
      return new DecompressingBinInputStream(*stream, true);
    }
    if (magic == "BZh")
    {
      return new DecompressingBinInputStream(*stream, false);
    }
    return new ArrowBinInputStream(*stream);
#else
    return nullptr; // not constructible without S3 support
#endif
  }

  bool S3InputSource::isS3URI(const std::string& uri)
  {
    return StringUtils::hasPrefix(uri, "s3://");
  }

  bool S3InputSource::isAvailable()
  {
#ifdef OPENMS_WITH_ARROW_S3
    return true;
#else
    return false;
#endif
  }

  Int64 S3InputSource::objectSize(const std::string& uri)
  {
    const S3InputSource source(uri);
#ifdef OPENMS_WITH_ARROW_S3
    auto size = source.impl_->object->GetSize();
    if (!size.ok())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, uri, size.status().ToString());
    }
    return *size;
#else
    return 0; // the constructor threw
#endif
  }

  std::string S3InputSource::readRange(const std::string& uri, Int64 offset, Int64 length)
  {
    const S3InputSource source(uri);
#ifdef OPENMS_WITH_ARROW_S3
    auto bytes = source.impl_->object->ReadAt(offset, length);
    if (!bytes.ok())
    {
      OPENMS_LOG_ERROR << "Error reading " << uri << ": " << bytes.status().ToString() << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, uri);
    }
    return (*bytes)->ToString();
#else
    (void)offset;
    (void)length;
    return {}; // the constructor threw
#endif
  }

} // namespace OpenMS
