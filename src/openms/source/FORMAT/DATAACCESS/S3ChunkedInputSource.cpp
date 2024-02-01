// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/DATAACCESS/S3ChunkedInputSource.h>
#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/HeadObjectRequest.h>


namespace OpenMS {

    S3ChunkedInputSource::S3ChunkedInputSource(const std::string& s3uri)
    {
        initializeAwsSdk_();
        parseS3Uri_(s3uri);
    }

    xercesc::BinInputStream* S3ChunkedInputSource::makeStream() const {
        Aws::Client::ClientConfiguration clientConfig;

        Aws::S3::S3Client s3Client(clientConfig);

        Aws::S3::Model::HeadObjectRequest request;
        request.WithBucket(m_bucketName).WithKey(m_objectKey);

        auto outcome = s3Client.HeadObject(request);
        std::string encoding;
        if (outcome.IsSuccess()) {
            auto metadata = outcome.GetResult().GetMetadata();
            auto encodingIter = metadata.find("Content-Encoding");
            if (encodingIter != metadata.end()) {
                encoding = encodingIter->second;
            }
            else
            {
                // First request to get the first magic bytes
                Aws::S3::Model::GetObjectRequest firstRequest;
                firstRequest.WithBucket(m_bucketName).WithKey(m_objectKey);
                firstRequest.SetRange("bytes=0-2");

                auto firstOutcome = s3Client.GetObject(firstRequest);

                if (firstOutcome.IsSuccess()) {
                    Aws::IOStream& firstStream = firstOutcome.GetResult().GetBody();
                    unsigned char header[3];
                    firstStream.read(reinterpret_cast<char*>(header), sizeof(header));
                    if (header[0] == 0x1F && header[1] == 0x8B) {
                        // The file is gzip compressed
                        encoding = "gzip";
                    } else if (header[0] == 'B' && header[1] == 'Z' && header[2] == 'h') {
                        // The file is bzip2 compressed
                        encoding = "bzip2";
                    }
                } else {
                    OPENMS_LOG_ERROR << "Error: AWS SDK GetObject for magic bytes: " <<
                        firstOutcome.GetError().GetExceptionName() << " " <<
                        firstOutcome.GetError().GetMessage() << std::endl;
                    return nullptr;
                }
            }
        } else {
            OPENMS_LOG_ERROR << "Error: AWS SDK HeadObject: " <<
                outcome.GetError().GetExceptionName() << " " <<
                outcome.GetError().GetMessage() << std::endl;
            return nullptr;
        }

        Aws::S3::Model::GetObjectRequest getObjectRequest;
        getObjectRequest.WithBucket(m_bucketName).WithKey(m_objectKey);

        auto getObjectOutcome = new Aws::S3::Model::GetObjectOutcome(s3Client.GetObject(getObjectRequest));
        if (getObjectOutcome->IsSuccess()) {
            if (encoding == "gzip") {
                return new S3ChunkedGzipBinInputStream(s3Client, getObjectRequest);
            } else if (encoding == "bzip2") {
                return new S3ChunkedBzip2BinInputStream(s3Client, getObjectRequest);
            } else {
                // The object is not compressed, or it is compressed with a different format
                return new S3ChunkedBinInputStream(s3Client, getObjectRequest);
            }
        } else {
            OPENMS_LOG_ERROR << "Error: AWS SDK GetObject: " <<
                getObjectOutcome->GetError().GetExceptionName() << " " <<
                getObjectOutcome->GetError().GetMessage() << std::endl;
            return nullptr;
        }
    }

    void S3ChunkedInputSource::initializeAwsSdk_() {
        Aws::SDKOptions options;
        Aws::InitAPI(options);
    }

    void S3ChunkedInputSource::cleanupAwsSdk_() {
        Aws::ShutdownAPI(Aws::SDKOptions());
    }

    void S3ChunkedInputSource::parseS3Uri_(std::string s3Uri) {
        // Remove the "s3://" prefix if present
        if (s3Uri.compare(0, 5, "s3://") == 0) {
            s3Uri = s3Uri.substr(5);
        }

        // Find the first occurrence of '/' character
        size_t slashPos = s3Uri.find('/');
        if (slashPos == std::string::npos) {
            // Invalid S3 URI format
            // Handle error
            return;
        }

        m_bucketName = s3Uri.substr(0, slashPos);
        m_objectKey = s3Uri.substr(slashPos + 1);
    }

    S3ChunkedBinInputStream::S3ChunkedBinInputStream(const Aws::S3::S3Client& client, const Aws::S3::Model::GetObjectRequest& req, unsigned long chunkSize)
        : m_req(req), m_client(client), m_position(0), m_chunkSize(chunkSize)
    {
        // load first chunk
        m_req.SetRange("bytes=0-" + std::to_string(m_chunkSize - 1));
        Aws::S3::Model::GetObjectOutcome outcome = m_client.GetObject(m_req);
        if (outcome.IsSuccess()) {
            m_currentChunk = outcome.GetResultWithOwnership();
        } else {
            OPENMS_LOG_ERROR << "Error: AWS SDK GetObject: " << outcome.GetError().GetExceptionName() << ": " << outcome.GetError().GetMessage() << std::endl;
        }
        m_currentChunkEnd = m_chunkSize - 1;
        std::string contentRange = m_currentChunk.GetContentRange();
        // Process the ContentRange header...
        // Example: bytes 0-524287/2000000
        size_t slashPos = contentRange.find('/');
        if (slashPos != std::string::npos) {
            std::string totalSizeStr = contentRange.substr(slashPos + 1);
            m_totalSize = std::stoull(totalSizeStr);
        } else {
            OPENMS_LOG_ERROR << "Error: AWS SDK GetObject: Invalid ContentRange header" << std::endl;
        }
    }

    XMLFilePos S3ChunkedBinInputStream::curPos() const { return m_position; }

    XMLSize_t S3ChunkedBinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) {
        auto currentBodyStream = &m_currentChunk.GetBody();
        unsigned long totalBytesRead = 0;
        // while maxToRead is not reached
        while (totalBytesRead < maxToRead) {
            // try to read from current chunk
            if (currentBodyStream->good())
            {
                currentBodyStream->read(reinterpret_cast<char*>(toFill + totalBytesRead), maxToRead - totalBytesRead);
                totalBytesRead += currentBodyStream->gcount();
                m_position += currentBodyStream->gcount();
            }
            else if (currentBodyStream->eof())
            {
                // if the end of the current chunk was already behind the total size, loading a new chunk is not necessary.
                // we're done.
                if (m_currentChunkEnd >= m_totalSize)
                {
                    return static_cast<XMLSize_t>(totalBytesRead);
                }
                // we reached the end of the current chunk
                // load next chunk
                m_req.SetRange("bytes=" + std::to_string(m_currentChunkEnd + 1) + "-" + std::to_string(m_currentChunkEnd + m_chunkSize));
                m_currentChunkEnd = m_currentChunkEnd + m_chunkSize;
                Aws::S3::Model::GetObjectOutcome outcome = m_client.GetObject(m_req);
                if (outcome.IsSuccess()) {
                    m_currentChunk = outcome.GetResultWithOwnership();
                    currentBodyStream = &m_currentChunk.GetBody();
                } else {
                    OPENMS_LOG_ERROR << "Error: " << outcome.GetError().GetExceptionName() << ": " << outcome.GetError().GetMessage() << std::endl;
                }
            }
            else
            {
                OPENMS_LOG_ERROR << "Error: AWS SDK result stream: " << currentBodyStream->rdstate() << std::endl;
            }
        }
        return static_cast<XMLSize_t>(totalBytesRead);
    }

    const XMLCh* S3ChunkedBinInputStream::getContentType() const { return nullptr; }

    S3ChunkedGzipBinInputStream::S3ChunkedGzipBinInputStream(const Aws::S3::S3Client& client, const Aws::S3::Model::GetObjectRequest& req, unsigned long chunkSize)
        : m_req(req), m_client(client), m_position(0), m_chunkSize(chunkSize), m_decompressedBuffer(), m_zStream()
    {
        // Create a GZIP decompression stream
        m_zStream.zalloc = Z_NULL;
        m_zStream.zfree = Z_NULL;
        m_zStream.opaque = Z_NULL;
        m_zStream.avail_in = 0;
        m_zStream.next_in = Z_NULL;

        int windowBits = 15 + 16; // Default windowBits for gzip
        if (inflateInit2(&m_zStream, windowBits) != Z_OK) {
            // Handle error
            throw new std::runtime_error("Error occurred during initializing Gzip stream. Zlib error message: " + std::string(m_zStream.msg));
        }

        // load first chunk
        m_req.SetRange("bytes=0-" + std::to_string(m_chunkSize - 1));
        Aws::S3::Model::GetObjectOutcome outcome = m_client.GetObject(m_req);
        if (outcome.IsSuccess()) {
            m_currentChunk = outcome.GetResultWithOwnership();
        } else {
            throw new std::runtime_error("Error getting chunk. AWS SDK error message: " + outcome.GetError().GetMessage());
        }
        m_currentChunkEnd = m_chunkSize - 1;
        std::string contentRange = m_currentChunk.GetContentRange();
        // Process the ContentRange header...
        // Example: bytes 0-524287/2000000
        size_t slashPos = contentRange.find('/');
        if (slashPos != std::string::npos) {
            std::string totalSizeStr = contentRange.substr(slashPos + 1);
            m_totalSize = std::stoull(totalSizeStr);
        } else {
            throw new std::runtime_error("Error: AWS SDK GetObject: Invalid ContentRange header: " + contentRange);
        }
    }

    S3ChunkedGzipBinInputStream::~S3ChunkedGzipBinInputStream() {
        inflateEnd(&m_zStream);
    }

    XMLFilePos S3ChunkedGzipBinInputStream::curPos() const { return m_position; }

    XMLSize_t S3ChunkedGzipBinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) {
        auto currentBodyStream = &m_currentChunk.GetBody();
        unsigned long totalBytesRead = 0;
        // while maxToRead is not reached
        while (totalBytesRead < maxToRead) {
            unsigned long remainingBytes = maxToRead - totalBytesRead;
            // read from current chunk
            if (currentBodyStream->good())
            {
                do {
                    m_zStream.avail_out = remainingBytes;
                    m_zStream.next_out = reinterpret_cast<Bytef*>(toFill + totalBytesRead);

                    if (m_zStream.avail_in == 0)
                    { // no bytes from last call in the zlib buffer anymore.
                        // Read more from current chunk if available
                        if (currentBodyStream->good()) {
                            currentBodyStream->read(reinterpret_cast<char*>(m_decompressedBuffer), 1024);
                            m_zStream.avail_in = currentBodyStream->gcount();
                            m_zStream.next_in = reinterpret_cast<Bytef*>(m_decompressedBuffer);
                        }
                        else {
                            // Handle error or end of stream
                            break;
                        }
                    }

                    int ret = inflate(&m_zStream, Z_NO_FLUSH);
                    if (ret == Z_STREAM_END) {
                        // Reached end of compressed data
                        unsigned long bytesRead = remainingBytes - m_zStream.avail_out;
                        totalBytesRead += bytesRead;
                        m_position += bytesRead;
                        remainingBytes = m_zStream.avail_out;
                        break;
                    }
                    else if (ret != Z_OK) {
                        OPENMS_LOG_ERROR << "Error occurred during decompression. Zlib error code: " << ret << std::endl;
                        break;
                    }

                    unsigned long bytesRead = remainingBytes - m_zStream.avail_out;
                    totalBytesRead += bytesRead;
                    m_position += bytesRead;
                    remainingBytes = m_zStream.avail_out;  // Update remainingBytes after each call to inflate
                } while (remainingBytes > 0);  // Continue the loop while there is still space in toFill
            }
            else if (currentBodyStream->eof())
            {
                // if the end of the current chunk was already behind the total size, loading a new chunk is not necessary
                if (m_currentChunkEnd >= m_totalSize)
                {
                    // We've reached the end of the file, but there might still be data in the zlib buffer.
                    m_zStream.avail_out = remainingBytes;
                    m_zStream.next_out = reinterpret_cast<Bytef*>(toFill + totalBytesRead);

                    int ret;
                    do { // emptying remaining data in zlib buffer
                        ret = inflate(&m_zStream, Z_NO_FLUSH);
                        if (ret == Z_STREAM_END) {
                            // Reached end of compressed data
                            unsigned long bytesRead = remainingBytes - m_zStream.avail_out;
                            totalBytesRead += bytesRead;
                            m_position += bytesRead;
                            remainingBytes = m_zStream.avail_out;
                            break;
                        }
                        else if (ret == Z_OK) {
                            // toFill is full
                            unsigned long bytesRead = remainingBytes - m_zStream.avail_out;
                            totalBytesRead += bytesRead;
                            m_position += bytesRead;
                            remainingBytes = m_zStream.avail_out;
                            break;
                        }
                        else if (ret != Z_BUF_ERROR) {
                            OPENMS_LOG_ERROR << "Error occurred during decompression. Zlib error code: " << ret << std::endl;
                            break;
                        }
                    } while (true);
                    return static_cast<XMLSize_t>(totalBytesRead);
                }
                // we reached the end of the current chunk
                // load next chunk
                m_req.SetRange("bytes=" + std::to_string(m_currentChunkEnd + 1) + "-" + std::to_string(m_currentChunkEnd + m_chunkSize));
                m_currentChunkEnd = m_currentChunkEnd + m_chunkSize;
                Aws::S3::Model::GetObjectOutcome outcome = m_client.GetObject(m_req);
                if (outcome.IsSuccess()) {
                    m_currentChunk = outcome.GetResultWithOwnership();
                    currentBodyStream = &m_currentChunk.GetBody();
                } else {
                    OPENMS_LOG_ERROR << "Error loading chunk from S3: " << outcome.GetError().GetExceptionName() << ": " << outcome.GetError().GetMessage() << std::endl;
                }
            }
            else
            { // neither good nor eof -> error
                OPENMS_LOG_ERROR << "Error in AWS S3 SDK result stream: " << currentBodyStream->rdstate() << std::endl;
            }
        }
        return static_cast<XMLSize_t>(totalBytesRead);
    }

    const XMLCh* S3ChunkedGzipBinInputStream::getContentType() const {return nullptr;};

    S3ChunkedBzip2BinInputStream::S3ChunkedBzip2BinInputStream(const Aws::S3::S3Client& client, const Aws::S3::Model::GetObjectRequest& req, unsigned long chunkSize)
        : m_req(req), m_client(client), m_position(0), m_chunkSize(chunkSize), m_decompressedBuffer(), m_bzStream()
    {
        // Create a GZIP decompression stream
        m_bzStream.bzalloc = NULL;
        m_bzStream.bzfree = NULL;
        m_bzStream.opaque = NULL;
        m_bzStream.avail_in = 0;
        m_bzStream.next_in = NULL;

        if (BZ2_bzDecompressInit(&m_bzStream, 0, 0) != BZ_OK) {
            throw new std::runtime_error("Error occurred during bzip decompression");
        }

        // load first chunk
        m_req.SetRange("bytes=0-" + std::to_string(m_chunkSize - 1));
        Aws::S3::Model::GetObjectOutcome outcome = m_client.GetObject(m_req);
        if (outcome.IsSuccess()) {
            m_currentChunk = outcome.GetResultWithOwnership();
        } else {
            throw new std::runtime_error("Error getting chunk. AWS SDK error message: " + outcome.GetError().GetMessage());
        }
        m_currentChunkEnd = m_chunkSize - 1;
        std::string contentRange = m_currentChunk.GetContentRange();
        // Process the ContentRange header...
        // Example: bytes 0-524287/2000000
        size_t slashPos = contentRange.find('/');
        if (slashPos != std::string::npos) {
            std::string totalSizeStr = contentRange.substr(slashPos + 1);
            m_totalSize = std::stoull(totalSizeStr);
        } else {
            throw new std::runtime_error("Error: AWS SDK GetObject: Invalid ContentRange header: " + contentRange);
        }
    }

    S3ChunkedBzip2BinInputStream::~S3ChunkedBzip2BinInputStream() {
        BZ2_bzDecompressEnd(&m_bzStream);
    }

    XMLFilePos S3ChunkedBzip2BinInputStream::curPos() const { return m_position; }

    XMLSize_t S3ChunkedBzip2BinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) {
        auto currentBodyStream = &m_currentChunk.GetBody();
        unsigned long totalBytesRead = 0;
        // while maxToRead is not reached
        while (totalBytesRead < maxToRead) {
            unsigned long remainingBytes = maxToRead - totalBytesRead;
            // read from current chunk
            if (currentBodyStream->good())
            {
                do {
                    m_bzStream.avail_out = remainingBytes;
                    m_bzStream.next_out = reinterpret_cast<char*>(toFill + totalBytesRead);

                    if (m_bzStream.avail_in == 0)
                    { // no bytes from last call in the zlib buffer anymore.
                        // Read more from current chunk if available
                        if (currentBodyStream->good()) {
                            currentBodyStream->read(reinterpret_cast<char*>(m_decompressedBuffer), 1024);
                            m_bzStream.avail_in = currentBodyStream->gcount();
                            m_bzStream.next_in = reinterpret_cast<char*>(m_decompressedBuffer);
                        }
                        else {
                            // Handle error or end of stream
                            break;
                        }
                    }

                    int ret = BZ2_bzDecompress(&m_bzStream);
                    if (ret == BZ_STREAM_END) {
                        // Reached end of compressed data
                        unsigned long bytesRead = remainingBytes - m_bzStream.avail_out;
                        totalBytesRead += bytesRead;
                        m_position += bytesRead;
                        remainingBytes = m_bzStream.avail_out;
                        break;
                    }
                    else if (ret != BZ_OK) {
                        OPENMS_LOG_ERROR << "Error occurred during decompression. Bzip2 error code: " << ret << std::endl;
                        break;
                    }

                    unsigned long bytesRead = remainingBytes - m_bzStream.avail_out;
                    totalBytesRead += bytesRead;
                    m_position += bytesRead;
                    remainingBytes = m_bzStream.avail_out;  // Update remainingBytes after each call to inflate
                } while (remainingBytes > 0);  // Continue the loop while there is still space in toFill
            }
            else if (currentBodyStream->eof())
            {
                // if the end of the current chunk was already behind the total size, loading a new chunk is not necessary
                if (m_currentChunkEnd >= m_totalSize)
                {
                    // We've reached the end of the file, but there might still be data in the zlib buffer.
                    m_bzStream.avail_out = remainingBytes;
                    m_bzStream.next_out = reinterpret_cast<char*>(toFill + totalBytesRead);

                    int ret;
                    do { // emptying remaining data in zlib buffer
                        ret = BZ2_bzDecompress(&m_bzStream);
                        if (ret == BZ_STREAM_END) {
                            // Reached end of compressed data
                            unsigned long bytesRead = remainingBytes - m_bzStream.avail_out;
                            totalBytesRead += bytesRead;
                            m_position += bytesRead;
                            remainingBytes = m_bzStream.avail_out;
                            break;
                        }
                        else if (ret == BZ_OK) {
                            // toFill is full
                            unsigned long bytesRead = remainingBytes - m_bzStream.avail_out;
                            totalBytesRead += bytesRead;
                            m_position += bytesRead;
                            remainingBytes = m_bzStream.avail_out;
                            break;
                        }
                        else if (ret != BZ_OK && ret != BZ_STREAM_END) {
                            OPENMS_LOG_ERROR << "Error occurred during decompression. Bzib2 error code: " << ret << std::endl;
                            break;
                        }
                    } while (true);
                    return static_cast<XMLSize_t>(totalBytesRead);
                }
                // we reached the end of the current chunk
                // load next chunk
                m_req.SetRange("bytes=" + std::to_string(m_currentChunkEnd + 1) + "-" + std::to_string(m_currentChunkEnd + m_chunkSize));
                m_currentChunkEnd = m_currentChunkEnd + m_chunkSize;
                Aws::S3::Model::GetObjectOutcome outcome = m_client.GetObject(m_req);
                if (outcome.IsSuccess()) {
                    m_currentChunk = outcome.GetResultWithOwnership();
                    currentBodyStream = &m_currentChunk.GetBody();
                } else {
                    OPENMS_LOG_ERROR << "Error loading chunk from S3: " << outcome.GetError().GetExceptionName() << ": " << outcome.GetError().GetMessage() << std::endl;
                }
            }
            else
            { // neither good nor eof -> error
                OPENMS_LOG_ERROR << "Error in AWS S3 SDK result stream: " << currentBodyStream->rdstate() << std::endl;
            }
        }
        return static_cast<XMLSize_t>(totalBytesRead);
    }

    const XMLCh* S3ChunkedBzip2BinInputStream::getContentType() const {return nullptr;};

}; // namespace OpenMS
