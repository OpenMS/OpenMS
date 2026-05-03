// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/S3InputSource.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/HeadObjectRequest.h>

#include <stdexcept>

namespace OpenMS
{

  // ---------------------------------------------------------------------------
  // S3InputSource
  // ---------------------------------------------------------------------------

  S3InputSource::S3InputSource(const std::string& s3uri)
  {
    AwsSdkHelper::initializeAwsSdk();
    parseS3Uri_(s3uri);
  }

  xercesc::BinInputStream* S3InputSource::makeStream() const
  {
    Aws::Client::ClientConfiguration client_config;
    Aws::S3::S3Client s3_client(client_config);

    // Determine compression from metadata or magic bytes
    std::string encoding;

    // Try Content-Encoding metadata first via HeadObject
    Aws::S3::Model::HeadObjectRequest head_request;
    head_request.WithBucket(bucket_name_).WithKey(object_key_);

    auto head_outcome = s3_client.HeadObject(head_request);
    if (head_outcome.IsSuccess())
    {
      auto metadata = head_outcome.GetResult().GetMetadata();
      auto it = metadata.find("Content-Encoding");
      if (it != metadata.end())
      {
        encoding = it->second;
      }
      else
      {
        // Check magic bytes by fetching the first 3 bytes
        Aws::S3::Model::GetObjectRequest magic_request;
        magic_request.WithBucket(bucket_name_).WithKey(object_key_);
        magic_request.SetRange("bytes=0-2");

        auto magic_outcome = s3_client.GetObject(magic_request);
        if (magic_outcome.IsSuccess())
        {
          Aws::IOStream& magic_stream = magic_outcome.GetResult().GetBody();
          unsigned char header[3] = {};
          magic_stream.read(reinterpret_cast<char*>(header), sizeof(header));

          if (header[0] == 0x1F && header[1] == 0x8B)
          {
            encoding = "gzip";
          }
          else if (header[0] == 'B' && header[1] == 'Z' && header[2] == 'h')
          {
            encoding = "bzip2";
          }
        }
        else
        {
          OPENMS_LOG_ERROR << "Error: AWS SDK GetObject for magic bytes: "
                          << magic_outcome.GetError().GetExceptionName() << " "
                          << magic_outcome.GetError().GetMessage() << std::endl;
          return nullptr;
        }
      }
    }
    else
    {
      OPENMS_LOG_ERROR << "Error: AWS SDK HeadObject: "
                      << head_outcome.GetError().GetExceptionName() << " "
                      << head_outcome.GetError().GetMessage() << std::endl;
      return nullptr;
    }

    // Fetch the full object
    Aws::S3::Model::GetObjectRequest get_request;
    get_request.WithBucket(bucket_name_).WithKey(object_key_);

    auto get_outcome = std::make_shared<Aws::S3::Model::GetObjectOutcome>(
      s3_client.GetObject(get_request));

    if (!get_outcome->IsSuccess())
    {
      OPENMS_LOG_ERROR << "Error: AWS SDK GetObject: "
                      << get_outcome->GetError().GetExceptionName() << " "
                      << get_outcome->GetError().GetMessage() << std::endl;
      return nullptr;
    }

    // Return the appropriate stream based on compression
    if (encoding == "gzip")
    {
      return new S3GzipBinInputStream(get_outcome);
    }
    else if (encoding == "bzip2")
    {
      return new S3Bzip2BinInputStream(get_outcome);
    }
    else
    {
      return new S3BinInputStream(get_outcome);
    }
  }

  void S3InputSource::parseS3Uri_(const std::string& s3uri)
  {
    std::string uri = s3uri;
    if (uri.compare(0, 5, "s3://") == 0)
    {
      uri = uri.substr(5);
    }

    Size slash_pos = uri.find('/');
    if (slash_pos == std::string::npos)
    {
      throw std::runtime_error("Invalid S3 URI format: missing '/' after bucket name in '" + s3uri + "'");
    }

    bucket_name_ = uri.substr(0, slash_pos);
    object_key_ = uri.substr(slash_pos + 1);

    if (bucket_name_.empty() || object_key_.empty())
    {
      throw std::runtime_error("Invalid S3 URI: bucket name or object key is empty in '" + s3uri + "'");
    }
  }

  // ---------------------------------------------------------------------------
  // S3BinInputStream
  // ---------------------------------------------------------------------------

  S3BinInputStream::S3BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome)
    : outcome_(std::move(outcome)),
      stream_(&outcome_->GetResult().GetBody()),
      position_(0)
  {
  }

  XMLFilePos S3BinInputStream::curPos() const
  {
    return position_;
  }

  XMLSize_t S3BinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead)
  {
    if (stream_ && toFill && stream_->good())
    {
      stream_->read(reinterpret_cast<char*>(toFill), maxToRead);
      std::streamsize bytes_read = stream_->gcount();
      position_ += bytes_read;
      return static_cast<XMLSize_t>(bytes_read);
    }
    return 0;
  }

  const XMLCh* S3BinInputStream::getContentType() const
  {
    return nullptr;
  }

  // ---------------------------------------------------------------------------
  // S3GzipBinInputStream
  // ---------------------------------------------------------------------------

  S3GzipBinInputStream::S3GzipBinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome)
    : outcome_(std::move(outcome)),
      stream_(&outcome_->GetResult().GetBody()),
      z_stream_(),
      buffer_(),
      position_(0)
  {
    z_stream_.zalloc = Z_NULL;
    z_stream_.zfree = Z_NULL;
    z_stream_.opaque = Z_NULL;
    z_stream_.avail_in = 0;
    z_stream_.next_in = Z_NULL;

    // windowBits = 15 + 16 enables gzip decoding
    int ret = inflateInit2(&z_stream_, 15 + 16);
    if (ret != Z_OK)
    {
      throw std::runtime_error("Error initializing gzip stream: " +
        std::string(z_stream_.msg ? z_stream_.msg : "unknown zlib error"));
    }
  }

  S3GzipBinInputStream::~S3GzipBinInputStream()
  {
    inflateEnd(&z_stream_);
  }

  XMLFilePos S3GzipBinInputStream::curPos() const
  {
    return position_;
  }

  XMLSize_t S3GzipBinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead)
  {
    z_stream_.avail_out = static_cast<uInt>(maxToRead);
    z_stream_.next_out = reinterpret_cast<Bytef*>(toFill);

    while (z_stream_.avail_out > 0)
    {
      // Refill input buffer when empty
      if (z_stream_.avail_in == 0)
      {
        if (!stream_ || !stream_->good())
        {
          break;
        }
        stream_->read(reinterpret_cast<char*>(buffer_), sizeof(buffer_));
        std::streamsize bytes_read = stream_->gcount();
        if (bytes_read == 0)
        {
          break;
        }
        z_stream_.avail_in = static_cast<uInt>(bytes_read);
        z_stream_.next_in = buffer_;
      }

      int ret = inflate(&z_stream_, Z_NO_FLUSH);
      if (ret == Z_STREAM_END)
      {
        break;
      }
      else if (ret != Z_OK)
      {
        OPENMS_LOG_ERROR << "Error during gzip decompression: "
                        << (z_stream_.msg ? z_stream_.msg : "unknown") << std::endl;
        break;
      }
    }

    XMLSize_t decompressed = maxToRead - z_stream_.avail_out;
    position_ += decompressed;
    return decompressed;
  }

  const XMLCh* S3GzipBinInputStream::getContentType() const
  {
    return nullptr;
  }

  // ---------------------------------------------------------------------------
  // S3Bzip2BinInputStream
  // ---------------------------------------------------------------------------

  S3Bzip2BinInputStream::S3Bzip2BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome)
    : outcome_(std::move(outcome)),
      stream_(&outcome_->GetResult().GetBody()),
      bz_stream_(),
      buffer_(),
      position_(0)
  {
    bz_stream_.bzalloc = nullptr;
    bz_stream_.bzfree = nullptr;
    bz_stream_.opaque = nullptr;
    bz_stream_.avail_in = 0;
    bz_stream_.next_in = nullptr;

    if (BZ2_bzDecompressInit(&bz_stream_, 0, 0) != BZ_OK)
    {
      throw std::runtime_error("Error initializing bzip2 decompression");
    }
  }

  S3Bzip2BinInputStream::~S3Bzip2BinInputStream()
  {
    BZ2_bzDecompressEnd(&bz_stream_);
  }

  XMLFilePos S3Bzip2BinInputStream::curPos() const
  {
    return position_;
  }

  XMLSize_t S3Bzip2BinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead)
  {
    bz_stream_.avail_out = static_cast<unsigned int>(maxToRead);
    bz_stream_.next_out = reinterpret_cast<char*>(toFill);

    while (bz_stream_.avail_out > 0)
    {
      // Refill input buffer when empty
      if (bz_stream_.avail_in == 0)
      {
        if (!stream_ || !stream_->good())
        {
          break;
        }
        stream_->read(buffer_, sizeof(buffer_));
        std::streamsize bytes_read = stream_->gcount();
        if (bytes_read == 0)
        {
          break;
        }
        bz_stream_.avail_in = static_cast<unsigned int>(bytes_read);
        bz_stream_.next_in = buffer_;
      }

      int ret = BZ2_bzDecompress(&bz_stream_);
      if (ret == BZ_STREAM_END)
      {
        break;
      }
      else if (ret != BZ_OK)
      {
        OPENMS_LOG_ERROR << "Error during bzip2 decompression, code: " << ret << std::endl;
        break;
      }
    }

    XMLSize_t decompressed = maxToRead - bz_stream_.avail_out;
    position_ += decompressed;
    return decompressed;
  }

  const XMLCh* S3Bzip2BinInputStream::getContentType() const
  {
    return nullptr;
  }

} // namespace OpenMS
