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
#include <aws/core/utils/memory/stl/AWSStreamFwd.h>
#include <aws/s3/model/GetObjectResult.h>
#include <aws/core/client/ClientConfiguration.h>
#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <zlib.h>
#include <bzlib.h>

#include <string>


namespace OpenMS {
    class S3InputSource : public xercesc::InputSource {
    public:
        S3InputSource(const std::string& s3uri);
        xercesc::BinInputStream* makeStream() const override;
    private:
        void initializeAwsSdk_();

        void cleanupAwsSdk_();

        void parseS3Uri_(std::string s3Uri);

        std::string m_bucketName;
        std::string m_objectKey;
    };

    class S3BinInputStream : public xercesc::BinInputStream {
    public:
        explicit S3BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> getObjectOutcome);

        XMLFilePos curPos() const override;

        XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

        const XMLCh* getContentType() const override;

    private:
        std::shared_ptr<Aws::S3::Model::GetObjectOutcome> m_getObjectOutcome;
        Aws::IOStream* m_stream;
        XMLFilePos m_position;
    };

    class S3GzipBinInputStream : public xercesc::BinInputStream {
    public:
        S3GzipBinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> getObjectOutcome);

        ~S3GzipBinInputStream();

        XMLFilePos curPos() const override;

        XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead) override;

        const XMLCh* getContentType() const override;

    private:
        std::shared_ptr<Aws::S3::Model::GetObjectOutcome> m_getObjectOutcome;
        Aws::IOStream* m_stream;
        z_stream m_zStream;
        Bytef m_buffer[1024];  // Decompressing buffer for reading from the S3 stream
        XMLFilePos m_position;
    };

    class S3Bzip2BinInputStream : public xercesc::BinInputStream
    {
    public:
        S3Bzip2BinInputStream(std::shared_ptr<Aws::S3::Model::GetObjectOutcome> getObjectOutcome);
        ~S3Bzip2BinInputStream();

        XMLFilePos curPos() const;
        XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead);
        const XMLCh* getContentType() const;

    private:
        std::shared_ptr<Aws::S3::Model::GetObjectOutcome> m_getObjectOutcome;
        Aws::IOStream* m_stream;
        bz_stream m_bzStream;
        char m_buffer[1024];  // Decompressing buffer for reading from the stream
        XMLFilePos m_position;
    };
}
