// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once
#include <xercesc/sax/InputSource.hpp>
#include <xercesc/util/BinInputStream.hpp>
#include <aws/s3/model/GetObjectResult.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <zlib.h>
#include <bzlib.h>

#include <string>


namespace OpenMS {
    class S3ChunkedInputSource : public xercesc::InputSource {
    public:
        S3ChunkedInputSource(const std::string& s3uri);
        xercesc::BinInputStream* makeStream() const override;
    private:
        void initializeAwsSdk_();

        void cleanupAwsSdk_();

        void parseS3Uri_(std::string s3Uri);

        std::string m_bucketName;
        std::string m_objectKey;
    };

    class S3ChunkedBinInputStream : public xercesc::BinInputStream {
    public:
        explicit S3ChunkedBinInputStream(
            const Aws::S3::S3Client& client,
            const Aws::S3::Model::GetObjectRequest& req,
            unsigned long chunkSize = 1024 * 1024 * 100 // 100 MB should be a good default for most cases
        );

        XMLFilePos curPos() const override;

        // Note: Typical maxToRead is 48 KB (48 * 1024)
        XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

        const XMLCh* getContentType() const override;

    private:
        Aws::S3::Model::GetObjectRequest m_req; ///< Copy of the request (bucket + key + range) where the range will be updated
        Aws::S3::S3Client m_client; ///< Client to use for the request (with credentials etc.)
        XMLFilePos m_position; ///< Current position in the stream important for curPos()
        unsigned long m_chunkSize; ///< Size of each chunk to request from S3
        unsigned long m_currentChunkEnd; ///< End of the current chunk
        unsigned long m_totalSize; ///< Total size of the object to download/stream
        Aws::S3::Model::GetObjectResult m_currentChunk; ///< Current chunk of the object as returned by S3
    };

    class S3ChunkedGzipBinInputStream : public xercesc::BinInputStream {
    public:
        S3ChunkedGzipBinInputStream(
            const Aws::S3::S3Client& client,
            const Aws::S3::Model::GetObjectRequest& req,
            unsigned long chunkSize = 1024 * 1024 * 100
        );

        ~S3ChunkedGzipBinInputStream();

        XMLFilePos curPos() const override;

        XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

        const XMLCh* getContentType() const override;

    private:
        Aws::S3::Model::GetObjectRequest m_req; ///< Copy of the request (bucket + key + range) where the range will be updated
        Aws::S3::S3Client m_client; ///< Client to use for the request (with credentials etc.)
        XMLFilePos m_position; ///< Current position in the stream important for curPos()
        unsigned long m_chunkSize; ///< Size of each chunk to request from S3
        unsigned long m_currentChunkEnd; ///< End of the current chunk
        unsigned long m_totalSize; ///< Total size of the object to download/stream
        Aws::S3::Model::GetObjectResult m_currentChunk; ///< Current chunk of the object as returned by S3
        XMLByte m_decompressedBuffer[1024]; ///< Current decompressed buffer
        z_stream m_zStream; ///< zlib stream for decompression
    };

    class S3ChunkedBzip2BinInputStream : public xercesc::BinInputStream {
    public:
        S3ChunkedBzip2BinInputStream(
            const Aws::S3::S3Client& client,
            const Aws::S3::Model::GetObjectRequest& req,
            unsigned long chunkSize = 1024 * 1024 * 100
        );

        ~S3ChunkedBzip2BinInputStream();

        XMLFilePos curPos() const override;

        XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

        const XMLCh* getContentType() const override;

    private:
        Aws::S3::Model::GetObjectRequest m_req; ///< Copy of the request (bucket + key + range) where the range will be updated
        Aws::S3::S3Client m_client; ///< Client to use for the request (with credentials etc.)
        XMLFilePos m_position; ///< Current position in the stream important for curPos()
        unsigned long m_chunkSize; ///< Size of each chunk to request from S3
        unsigned long m_currentChunkEnd; ///< End of the current chunk
        unsigned long m_totalSize; ///< Total size of the object to download/stream
        Aws::S3::Model::GetObjectResult m_currentChunk; ///< Current chunk of the object as returned by S3
        XMLByte m_decompressedBuffer[1024]; ///< Current decompressed buffer
        bz_stream m_bzStream; ///< zlib stream for decompression
    };
}
