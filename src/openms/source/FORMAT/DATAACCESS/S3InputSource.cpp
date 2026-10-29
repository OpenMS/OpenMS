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
#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h>
#include <OpenMS/FORMAT/DecompressingInputStream.h>

#ifdef OPENMS_WITH_ARROW_S3
#include <arrow/buffer.h>
#include <arrow/filesystem/s3fs.h>
#include <arrow/io/buffered.h>
#include <arrow/io/interfaces.h>
#include <arrow/memory_pool.h>
#endif

#include <cstdlib>
#include <mutex>
#include <string_view>

namespace OpenMS
{

  struct S3InputSource::Impl
  {
    std::string uri;
    std::string name; ///< the URI for messages, without credentials (HTTPInputSource::redactURL())
#ifdef OPENMS_WITH_ARROW_S3
    std::shared_ptr<arrow::io::RandomAccessFile> object;
#endif
  };

  namespace
  {
#ifdef OPENMS_WITH_ARROW_S3
    // Each read of an S3 object is a request: the streams read the object in blocks of this size
    constexpr int64_t READ_BLOCK_SIZE = 8 * 1024 * 1024;

    // Arrow's S3 APIs are initialized before their first use and finalized at exit. Without the finalization, the
    // process can crash at exit (e.g. on macOS), also when S3 was initialized but no object was opened (an invalid URI).
    // The finalization uses statics that Arrow creates with its first S3 client, and exit handlers and static
    // destructors run in the reverse order of their registration. So a first client is made right away (it does not
    // connect), and the exit handler is registered after it, to run before those statics are destroyed.
    void ensureS3Initialized()
    {
      static std::once_flag once;
      static arrow::Status status;
      std::call_once(once, [] {
        status = arrow::fs::EnsureS3Initialized();
        if (!status.ok()) { return; }
        status = arrow::fs::S3FileSystem::Make(arrow::fs::S3Options::Anonymous()).status();
        std::atexit([] { (void)arrow::fs::EnsureS3Finalized(); });
      });
      if (!status.ok())
      {
        OPENMS_LOG_ERROR << "Cannot initialize S3: " << HTTPInputSource::redactURLs(status.ToString()) << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot initialize S3: " + HTTPInputSource::redactURLs(status.ToString()));
      }
    }

    // Opens the S3 object @p uri
    std::shared_ptr<arrow::io::RandomAccessFile> openObject(const std::string& uri)
    {
      ensureS3Initialized();
      const std::string name = HTTPInputSource::redactURL(uri); // the URI may hold an access key and secret
      std::string path;
      auto options = arrow::fs::S3Options::FromUri(uri, &path);
      // the errors are logged, as callers may rethrow them as other exceptions (e.g. MzMLFile::load(): a ParseError)
      if (!options.ok())
      {
        OPENMS_LOG_ERROR << "Invalid S3 URI " << name << ": " << HTTPInputSource::redactURLs(options.status().ToString(), uri) << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Invalid S3 URI '" + name + "': " + HTTPInputSource::redactURLs(options.status().ToString(), uri));
      }
      auto filesystem = arrow::fs::S3FileSystem::Make(*options);
      if (!filesystem.ok())
      {
        OPENMS_LOG_ERROR << "Cannot open " << name << ": " << HTTPInputSource::redactURLs(filesystem.status().ToString(), uri) << std::endl;
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name, HTTPInputSource::redactURLs(filesystem.status().ToString(), uri));
      }
      auto object = (*filesystem)->OpenInputFile(path);
      if (!object.ok())
      {
        OPENMS_LOG_ERROR << "Cannot open " << name << ": " << HTTPInputSource::redactURLs(object.status().ToString(), uri) << std::endl;
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name, HTTPInputSource::redactURLs(object.status().ToString(), uri));
      }
      return *object;
    }

    // The bytes of an Arrow stream, for the parser
    class ArrowInputStream :
      public InputStream
    {
public:
      ArrowInputStream(std::shared_ptr<arrow::io::InputStream> stream, const std::string& uri) :
        stream_(std::move(stream)),
        uri_(uri),
        name_(HTTPInputSource::redactURL(uri))
      {
      }

      UInt64 curPos() const override
      {
        return position_;
      }

      Size readBytes(unsigned char* to_fill, Size max_to_read) override
      {
        auto bytes_read = stream_->Read(static_cast<int64_t>(max_to_read), to_fill);
        if (!bytes_read.ok())
        {
          OPENMS_LOG_ERROR << "Error reading " << name_ << ": " << HTTPInputSource::redactURLs(bytes_read.status().ToString(), uri_) << std::endl;
          throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name_);
        }
        position_ += static_cast<UInt64>(*bytes_read);
        return static_cast<Size>(*bytes_read);
      }

private:
      std::shared_ptr<arrow::io::InputStream> stream_;
      std::string uri_;
      std::string name_; ///< the URI for messages, without credentials
      UInt64 position_ = 0;
    };
#endif
  } // namespace

  S3InputSource::S3InputSource(const std::string& uri) :
    impl_(std::make_unique<Impl>())
  {
    if (!isS3URI(uri))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Not an S3 URI (s3://bucket/key): '" + HTTPInputSource::redactURL(uri) + "'");
    }
    if (!isAvailable())
    {
      const std::string message = "Cannot read '" + HTTPInputSource::redactURL(uri) + "': OpenMS was built without S3 support (its Arrow library has none, "
                                  "e.g. vcpkg's arrow without the feature 's3').";
      OPENMS_LOG_ERROR << message << std::endl; // callers may rethrow it as another exception (e.g. MzMLFile::load())
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
    }
    impl_->uri = uri;
    impl_->name = HTTPInputSource::redactURL(uri);
#ifdef OPENMS_WITH_ARROW_S3
    impl_->object = openObject(uri);
#endif
  }

  S3InputSource::~S3InputSource() = default;

  std::unique_ptr<InputStream> S3InputSource::makeStream() const
  {
#ifdef OPENMS_WITH_ARROW_S3
    auto fail = [this](const arrow::Status& status) -> std::unique_ptr<InputStream>
    {
      OPENMS_LOG_ERROR << "Error reading " << impl_->name << ": " << HTTPInputSource::redactURLs(status.ToString(), impl_->uri) << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, impl_->name);
    };
    // a stream of its own for each call, read in large blocks
    auto size = impl_->object->GetSize();
    if (!size.ok()) { return fail(size.status()); }
    auto raw = arrow::io::RandomAccessFile::GetStream(impl_->object, 0, *size);
    if (!raw.ok()) { return fail(raw.status()); }
    auto stream = arrow::io::BufferedInputStream::Create(READ_BLOCK_SIZE, arrow::default_memory_pool(), *raw);
    if (!stream.ok()) { return fail(stream.status()); }
    // gzip or bzip2 compressed: by the magic bytes
    auto head = (*stream)->Peek(3);
    const std::string_view magic = head.ok() ? *head : std::string_view();
    const auto format = DecompressingInputStream::detect(reinterpret_cast<const unsigned char*>(magic.data()), magic.size());
    auto plain = std::make_unique<ArrowInputStream>(*stream, impl_->uri);
    if (format)
    {
      return std::make_unique<DecompressingInputStream>(std::move(plain), *format, impl_->name);
    }
    return plain;
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
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, source.impl_->name, HTTPInputSource::redactURLs(size.status().ToString(), uri));
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
      OPENMS_LOG_ERROR << "Error reading " << source.impl_->name << ": " << HTTPInputSource::redactURLs(bytes.status().ToString(), uri) << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, source.impl_->name);
    }
    return (*bytes)->ToString();
#else
    (void)offset;
    (void)length;
    return {}; // the constructor threw
#endif
  }

} // namespace OpenMS
