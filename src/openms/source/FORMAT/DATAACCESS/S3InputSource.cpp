// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/S3InputSource.h>
#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/HeadObjectRequest.h>

namespace OpenMS {

    S3InputSource::S3InputSource(const std::string& s3uri)
    {
        initializeAwsSdk_();
        parseS3Uri_(s3uri);
    }

    xercesc::BinInputStream* S3InputSource::makeStream() const {
        Aws::Client::ClientConfiguration clientConfig;
        //clientConfig.region = Aws::Region::US_EAST_1;

        // Create the default AWS SDK client
        //Aws::Auth::DefaultAWSCredentialsProviderChain credentialsProvider;

        // If the credentials and region are set in the AWS config file, they will be loaded here
        //clientConfig.credentialsProvider = Aws::MakeShared<Aws::Auth::DefaultAWSCredentialsProviderChain>("AllocationTag", credentialsProvider);

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
                    // Handle error here
                    std::cerr << "Error: AWS SDK GetObject for magic bytes: " <<
                        firstOutcome.GetError().GetExceptionName() << " " <<
                        firstOutcome.GetError().GetMessage() << std::endl;
                    return nullptr;
                }
            }
        } else {
            // Handle error here
            std::cerr << "Error: AWS SDK HeadObject: " <<
                outcome.GetError().GetExceptionName() << " " <<
                outcome.GetError().GetMessage() << std::endl;
            return nullptr;
        }

        Aws::S3::Model::GetObjectRequest getObjectRequest;
        getObjectRequest.WithBucket(m_bucketName).WithKey(m_objectKey);

        auto getObjectOutcome = new Aws::S3::Model::GetObjectOutcome(s3Client.GetObject(getObjectRequest));
        if (getObjectOutcome->IsSuccess()) {
            if (encoding == "gzip") {
                // The object is gzipped
                std::cout << "The file is gzip compressed" << std::endl;
                return new S3GzipBinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome>(getObjectOutcome));
            } else if (encoding == "bzip2") {
                // The object is bzipped
                std::cout << "The file is bzip compressed" << std::endl;
                return new S3Bzip2BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome>(getObjectOutcome));
            } else {
                // The object is not compressed, or it is compressed with a different format
                std::cout << "The file is uncompressed" << std::endl;
                return new S3BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome>(getObjectOutcome));
            }
            /* Reading magic bytes does not work since it consumes the stream -.-
            // Read the first few bytes of the file
            unsigned char buffer[3];
            getObjectOutcome->GetResult().GetBody().read(reinterpret_cast<char*>(buffer), sizeof(buffer));
            if (buffer[0] == 0x1F && buffer[1] == 0x8B) {
                // The file is gzip compressed
                std::cout << "The file is gzip compressed" << std::endl;
                return new S3GzipBinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome>(getObjectOutcome));
            } else if (buffer[0] == 'B' && buffer[1] == 'Z' && buffer[2] == 'h') {
                // The file is bzip2 compressed
                std::cout << "The file is bzip2 compressed" << std::endl;
                return new S3Bzip2BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome>(getObjectOutcome));
            } else {
                std::cout << "The file is uncompressed" << std::endl;
                return new S3BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome>(getObjectOutcome));
            }
            */
        } else {
            // Handle error here
            std::cerr << "Error: AWS SDK GetObject: " <<
                getObjectOutcome->GetError().GetExceptionName() << " " <<
                getObjectOutcome->GetError().GetMessage() << std::endl;
            return nullptr;
        }
    }

    void S3InputSource::initializeAwsSdk_() {
        Aws::SDKOptions options;
        Aws::InitAPI(options);
    }

    void S3InputSource::cleanupAwsSdk_() {
        Aws::ShutdownAPI(Aws::SDKOptions());
    }

    void S3InputSource::parseS3Uri_(std::string s3Uri) {
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

    /*
    void S3InputSource::readAwsConfiguration_() {
        m_clientConfig.region = Aws::Environment::GetEnv("AWS_REGION");
        m_clientConfig.endpointOverride = Aws::Environment::GetEnv("AWS_ENDPOINT");
        m_clientConfig.scheme = Aws::Environment::GetEnv("AWS_SCHEME");
        m_clientConfig.verifySSL = Aws::Utils::StringUtils::ConvertToBool(Aws::Environment::GetEnv("AWS_VERIFY_SSL"));
        // Set other configuration options as needed
    }*/

    S3BinInputStream::S3BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> getObjectOutcome)
    : m_getObjectOutcome(getObjectOutcome), m_stream(&getObjectOutcome->GetResult().GetBody()), m_position(0) {}

    XMLFilePos S3BinInputStream::curPos() const { return m_position; }

    XMLSize_t S3BinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) {

        if (m_stream && toFill) {
            if (m_stream->good()) {
                m_stream->read(reinterpret_cast<char*>(toFill), maxToRead);
                std::streamsize bytesRead = m_stream->gcount();
                m_position += bytesRead;
                return static_cast<XMLSize_t>(bytesRead);
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

    const XMLCh* S3BinInputStream::getContentType() const { return nullptr; }



    S3GzipBinInputStream::S3GzipBinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> getObjectOutcome)
        : m_getObjectOutcome(getObjectOutcome), m_stream(&getObjectOutcome->GetResult().GetBody()), m_zStream(), m_position(0)
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
            throw new std::runtime_error("Error occurred during decompression. Zlib error message: " + std::string(m_zStream.msg));
        }
        
    }

    S3GzipBinInputStream::~S3GzipBinInputStream() {
        inflateEnd(&m_zStream);
    }

    XMLFilePos S3GzipBinInputStream::curPos() const { return m_position; }

    XMLSize_t S3GzipBinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) {
            m_zStream.avail_out = static_cast<uInt>(maxToRead);
            m_zStream.next_out = reinterpret_cast<Bytef*>(toFill);

            // although maxToRead is usually 48kB we cannot know for sure
            // and must loop over smaller fixed buffer sizes for decompression
            while (m_zStream.avail_out > 0) {
                // Read more data from the stream
                if (m_zStream.avail_in == 0) {
                    m_stream->read(reinterpret_cast<char*>(m_buffer), sizeof(m_buffer));
                    std::streamsize bytesRead = m_stream->gcount();
                    m_position += bytesRead;

                    m_zStream.avail_in = static_cast<uInt>(bytesRead);
                    m_zStream.next_in = m_buffer;
                }

                int result = inflate(&m_zStream, Z_NO_FLUSH);
                if (result == Z_STREAM_END) {
                    // Reached the end of the compressed data
                    break;
                } else if (result != Z_OK) {
                    // Error occurred during decompression
                    // Handle error
                    std::cerr << "Error occurred during decompression. Zlib error message: " + std::string(m_zStream.msg) << std::endl;
                    return 0;
                }
            }

            XMLSize_t bytesRead = maxToRead - m_zStream.avail_out;
            return bytesRead;
        }

    const XMLCh* S3GzipBinInputStream::getContentType() const {return nullptr;};

    S3Bzip2BinInputStream::S3Bzip2BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> getObjectOutcome)
        : m_getObjectOutcome(getObjectOutcome), m_stream(&getObjectOutcome->GetResult().GetBody()), m_bzStream(), m_position(0)
    {
        m_bzStream.bzalloc = NULL;
        m_bzStream.bzfree = NULL;
        m_bzStream.opaque = NULL;
        m_bzStream.avail_in = 0;
        m_bzStream.next_in = NULL;

        if (BZ2_bzDecompressInit(&m_bzStream, 0, 0) != BZ_OK) {
            throw new std::runtime_error("Error occurred during decompression");
        }
    }

    S3Bzip2BinInputStream::~S3Bzip2BinInputStream() {
        BZ2_bzDecompressEnd(&m_bzStream);
    }

    XMLFilePos S3Bzip2BinInputStream::curPos() const { return m_position; }

    XMLSize_t S3Bzip2BinInputStream::readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) {
        m_bzStream.avail_out = static_cast<unsigned int>(maxToRead);
        m_bzStream.next_out = reinterpret_cast<char*>(toFill);

        while (m_bzStream.avail_out > 0) {
            if (m_bzStream.avail_in == 0) {
                m_stream->read(m_buffer, sizeof(m_buffer));
                std::streamsize bytesRead = m_stream->gcount();
                m_position += bytesRead;

                m_bzStream.avail_in = static_cast<unsigned int>(bytesRead);
                m_bzStream.next_in = m_buffer;
            }

            int result = BZ2_bzDecompress(&m_bzStream);
            if (result == BZ_STREAM_END) {
                break;
            } else if (result != BZ_OK) {
                std::cerr << "Error occurred during decompression" << std::endl;
                return 0;
            }
        }

        XMLSize_t bytesRead = maxToRead - m_bzStream.avail_out;
        return bytesRead;
    }

    const XMLCh* S3Bzip2BinInputStream::getContentType() const {return nullptr;}

}; // namespace OpenMS
