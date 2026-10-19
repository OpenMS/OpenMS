// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "S3InputSource.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileNameUtils.h>

#include <arrow/filesystem/filesystem.h>
#include <arrow/io/compressed.h>
#include <arrow/io/interfaces.h>
#include <arrow/buffer.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <arrow/util/compression.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>

namespace OpenMS
{
  namespace
  {
    template <typename T>
    T checked(arrow::Result<T> result, const char* operation)
    {
      if (!result.ok())
      {
        // Arrow messages may echo credentials or signed URI query parameters.
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    std::string("Arrow remote XML ") + operation + " failed ("
                                    + result.status().CodeAsString() + ")");
      }
      return std::move(result).ValueOrDie();
    }

    // Replay the magic bytes without seeking or reopening a remote object.
    class PrefixedStream final : public arrow::io::InputStream
    {
    public:
      explicit PrefixedStream(std::shared_ptr<arrow::io::InputStream> input) : input_(std::move(input))
      {
        while (size_ < 3)
        {
          const auto n = checked(input_->Read(3 - size_, prefix_.data() + size_), "compression-header read");
          if (n == 0) break;
          size_ += n;
        }
      }

      arrow::Status Close() override { return input_->Close(); }
      bool closed() const override { return input_->closed(); }
      arrow::Result<int64_t> Tell() const override { return position_; }
      arrow::Result<int64_t> Read(int64_t nbytes, void* out) override
      {
        if (nbytes < 0) return arrow::Status::Invalid("Negative read size");
        if (closed()) return arrow::Status::Invalid("Read from closed stream");
        const auto n = std::min(nbytes, size_ - offset_);
        if (n > 0)
        {
          std::memcpy(out, prefix_.data() + offset_, static_cast<size_t>(n));
          offset_ += n;
          position_ += n;
          return n;
        }
        auto result = input_->Read(nbytes, out);
        if (result.ok()) position_ += *result;
        return result;
      }
      arrow::Result<std::shared_ptr<arrow::Buffer>> Read(int64_t nbytes) override
      {
        auto buffer = arrow::AllocateResizableBuffer(nbytes);
        if (!buffer.ok()) return buffer.status();
        auto n = Read(nbytes, (*buffer)->mutable_data());
        if (!n.ok()) return n.status();
        auto status = (*buffer)->Resize(*n);
        if (!status.ok()) return status;
        return std::shared_ptr<arrow::Buffer>(std::move(*buffer));
      }
      arrow::Compression::type compression() const
      {
        if (size_ >= 2 && prefix_[0] == 0x1f && prefix_[1] == 0x8b) return arrow::Compression::GZIP;
        if (size_ == 3 && prefix_[0] == 'B' && prefix_[1] == 'Z' && prefix_[2] == 'h') return arrow::Compression::BZ2;
        return arrow::Compression::UNCOMPRESSED;
      }

    private:
      std::shared_ptr<arrow::io::InputStream> input_;
      std::array<unsigned char, 3> prefix_{};
      int64_t size_ = 0;
      int64_t offset_ = 0;
      int64_t position_ = 0;
    };
  }

  struct ArrowInputStream::Impl
  {
    std::shared_ptr<arrow::fs::FileSystem> filesystem;
    std::shared_ptr<arrow::io::InputStream> raw;
    std::unique_ptr<arrow::util::Codec> codec;
    std::shared_ptr<arrow::io::InputStream> stream;
    XMLFilePos position = 0;

    ~Impl()
    {
      if (stream) (void)stream->Close();
      if (raw && !raw->closed()) (void)raw->Close();
    }
  };

  ArrowInputStream::ArrowInputStream(std::shared_ptr<arrow::io::InputStream> stream,
                                   std::shared_ptr<arrow::fs::FileSystem> filesystem)
    : impl_(std::make_unique<Impl>())
  {
    if (!stream)
    {
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Null Arrow input stream");
    }
    impl_->filesystem = std::move(filesystem);
    impl_->raw = std::move(stream);
    auto prefix = std::make_shared<PrefixedStream>(impl_->raw);
    const auto compression = prefix->compression();
    if (compression == arrow::Compression::UNCOMPRESSED)
    {
      impl_->stream = std::move(prefix);
    }
    else
    {
      impl_->codec = checked(arrow::util::Codec::Create(compression), "compression codec initialization");
      impl_->stream = checked(arrow::io::CompressedInputStream::Make(impl_->codec.get(), prefix), "decompression setup");
    }
  }

  ArrowInputStream::~ArrowInputStream() = default;

  XMLFilePos ArrowInputStream::curPos() const { return impl_->position; }

  XMLSize_t ArrowInputStream::readBytes(XMLByte* const toFill, XMLSize_t maxToRead)
  {
    if (maxToRead == 0) return 0;
    if (!toFill)
    {
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Null XML destination buffer");
    }
    const auto n = checked(impl_->stream->Read(std::min<XMLSize_t>(maxToRead, 65536), toFill), "read/decompression");
    if (n < 0 || static_cast<uint64_t>(n) > std::min<XMLSize_t>(maxToRead, 65536)
        || static_cast<uint64_t>(n) > std::numeric_limits<XMLFilePos>::max() - impl_->position)
    {
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid Arrow XML read position");
    }
    impl_->position += n;
    return static_cast<XMLSize_t>(n);
  }

  const XMLCh* ArrowInputStream::getContentType() const { return nullptr; }

  S3InputSource::S3InputSource(const std::string& uri)
    : xercesc::InputSource(FileNameUtils::redactRemoteResource(uri).c_str()), uri_(uri)
  {
    if (uri.compare(0, 5, "s3://") != 0)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "", "Expected an s3:// URI");
    }
  }

  xercesc::BinInputStream* S3InputSource::makeStream() const
  {
    std::string path;
    auto filesystem = checked(arrow::fs::FileSystemFromUri(uri_, &path), "URI resolution (check URI options and S3 support)");
    auto input = checked(filesystem->OpenInputStream(path), "object open (check connectivity, credentials and permissions)");
    return new ArrowInputStream(std::move(input), std::move(filesystem));
  }
}
