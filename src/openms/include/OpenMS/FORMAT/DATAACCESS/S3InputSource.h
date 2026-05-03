// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Julianus Pfeuffer, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#ifdef WITH_S3

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/config.h>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/util/BinInputStream.hpp>

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectResult.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/core/utils/memory/stl/AWSStreamFwd.h>
#include <zlib.h>
#include <bzlib.h>

#include <memory>
#include <mutex>
#include <string>

namespace OpenMS
{

  /// @brief Helper to ensure the AWS SDK is initialized exactly once across all translation units
  namespace AwsSdkHelper
  {
    /// @brief Thread-safe one-time initialization of the AWS SDK
    inline void initializeAwsSdk()
    {
      static std::once_flag flag;
      static Aws::SDKOptions options;
      std::call_once(flag, []() {
        Aws::InitAPI(options);
        std::atexit([]() { Aws::ShutdownAPI(options); });
      });
    }
  } // namespace AwsSdkHelper

  /**
    @brief Xerces InputSource for reading XML files from Amazon S3.

    Provides streaming access to S3 objects for the Xerces XML parser.
    Automatically detects and handles gzip and bzip2 compressed objects.

    The S3 URI format is: s3://bucket-name/object-key

    Requires that the AWS SDK is available and configured (environment variables
    AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_DEFAULT_REGION, etc.).

    @ingroup Format
  */
  class OPENMS_DLLAPI S3InputSource :
    public xercesc::InputSource
  {
public:
    /**
      @brief Construct from an S3 URI.

      @param[in] s3uri S3 URI in the form s3://bucket/key

      @throws std::runtime_error if the URI format is invalid
    */
    explicit S3InputSource(const std::string& s3uri);

    /**
      @brief Create a BinInputStream for streaming data from S3.

      Detects compression by checking Content-Encoding metadata or magic bytes,
      and returns the appropriate decompressing stream.

      @return Pointer to a new BinInputStream (caller takes ownership), or nullptr on error
    */
    xercesc::BinInputStream* makeStream() const override;

private:
    /// Parse an s3:// URI into bucket and key components
    void parseS3Uri_(const std::string& s3uri);

    std::string bucket_name_; ///< S3 bucket name
    std::string object_key_;  ///< S3 object key
  };

  /**
    @brief BinInputStream that reads an S3 object as a byte stream.

    Wraps an AWS S3 GetObjectResult and provides sequential read access
    suitable for XML parsing.

    @ingroup Format
  */
  class OPENMS_DLLAPI S3BinInputStream :
    public xercesc::BinInputStream
  {
public:
    /**
      @brief Construct from a completed S3 GetObject outcome.

      @param[in] outcome Shared pointer to a successful GetObjectOutcome
    */
    explicit S3BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome);

    /// @brief Return the current byte position in the stream
    XMLFilePos curPos() const override;

    /**
      @brief Read up to @p maxToRead bytes into @p toFill.

      @param[out] toFill   Buffer to write data into
      @param[in]  maxToRead Maximum number of bytes to read

      @return Number of bytes actually read (0 at end-of-stream)
    */
    XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

    /// @brief Return the content type (always returns nullptr)
    const XMLCh* getContentType() const override;

private:
    std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome_; ///< Owns the GetObject result
    Aws::IOStream* stream_;  ///< Pointer to the result body stream (non-owning)
    XMLFilePos position_;    ///< Current read position
  };

  /**
    @brief BinInputStream that reads and decompresses gzip-encoded S3 objects.

    @ingroup Format
  */
  class OPENMS_DLLAPI S3GzipBinInputStream :
    public xercesc::BinInputStream
  {
public:
    /**
      @brief Construct from a completed S3 GetObject outcome containing gzip data.

      @param[in] outcome Shared pointer to a successful GetObjectOutcome

      @throws std::runtime_error if zlib initialization fails
    */
    explicit S3GzipBinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome);

    ~S3GzipBinInputStream() override;

    /// @brief Return the current byte position in the decompressed stream
    XMLFilePos curPos() const override;

    /**
      @brief Read and decompress up to @p maxToRead bytes into @p toFill.

      @param[out] toFill   Buffer to write decompressed data into
      @param[in]  maxToRead Maximum number of bytes to read

      @return Number of decompressed bytes written (0 at end-of-stream)
    */
    XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

    /// @brief Return the content type (always returns nullptr)
    const XMLCh* getContentType() const override;

private:
    std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome_; ///< Owns the GetObject result
    Aws::IOStream* stream_;      ///< Pointer to the result body stream (non-owning)
    z_stream z_stream_;          ///< zlib decompression state
    Bytef buffer_[4096];         ///< Compressed data read buffer
    XMLFilePos position_;        ///< Current decompressed byte position
  };

  /**
    @brief BinInputStream that reads and decompresses bzip2-encoded S3 objects.

    @ingroup Format
  */
  class OPENMS_DLLAPI S3Bzip2BinInputStream :
    public xercesc::BinInputStream
  {
public:
    /**
      @brief Construct from a completed S3 GetObject outcome containing bzip2 data.

      @param[in] outcome Shared pointer to a successful GetObjectOutcome

      @throws std::runtime_error if bzip2 initialization fails
    */
    explicit S3Bzip2BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome);

    ~S3Bzip2BinInputStream() override;

    /// @brief Return the current byte position in the decompressed stream
    XMLFilePos curPos() const override;

    /**
      @brief Read and decompress up to @p maxToRead bytes into @p toFill.

      @param[out] toFill   Buffer to write decompressed data into
      @param[in]  maxToRead Maximum number of bytes to read

      @return Number of decompressed bytes written (0 at end-of-stream)
    */
    XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

    /// @brief Return the content type (always returns nullptr)
    const XMLCh* getContentType() const override;

private:
    std::shared_ptr<Aws::S3::Model::GetObjectOutcome> outcome_; ///< Owns the GetObject result
    Aws::IOStream* stream_;       ///< Pointer to the result body stream (non-owning)
    bz_stream bz_stream_;         ///< bzip2 decompression state
    char buffer_[4096];           ///< Compressed data read buffer
    XMLFilePos position_;         ///< Current decompressed byte position
  };

} // namespace OpenMS

#endif // WITH_S3
