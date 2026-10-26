// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/DATAACCESS/S3InputSource.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <cstdlib>

using namespace OpenMS;

START_TEST(S3InputSource, "$Id$")

START_SECTION((static bool isS3URI(const std::string& uri)))
{
  TEST_TRUE(S3InputSource::isS3URI("s3://bucket/data/run.mzML"))
  TEST_TRUE(S3InputSource::isS3URI("s3://bucket/run.mzML?region=eu-central-1"))
  TEST_FALSE(S3InputSource::isS3URI("https://bucket.s3.amazonaws.com/run.mzML"))
  TEST_FALSE(S3InputSource::isS3URI("data/s3://run.mzML"))
}
END_SECTION

START_SECTION((explicit S3InputSource(const std::string& uri)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, S3InputSource("https://example.org/run.mzML"))
  TEST_EXCEPTION(Exception::IllegalArgument, S3InputSource(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML")))
  if (!S3InputSource::isAvailable())
  {
    // built without S3 support: s3:// URIs are an error, not a missing file
    TEST_EXCEPTION(Exception::IllegalArgument, S3InputSource("s3://bucket/run.mzML"))
    TEST_EXCEPTION(Exception::IllegalArgument, S3InputSource::objectSize("s3://bucket/run.mzML"))
  }
}
END_SECTION

START_SECTION((static bool isAvailable()))
{
  NOT_TESTABLE // depends on the Arrow library OpenMS was built with
}
END_SECTION

// Reading from S3 needs a server. Set OPENMS_TEST_S3_URI to the S3 URI of a copy of MzMLFile_1.mzML, e.g.
// s3://bucket/MzMLFile_1.mzML?endpoint_override=localhost:9000&scheme=http for a local S3-compatible server (with its
// credentials in AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY), to compare what is read from S3 with the local file.
START_SECTION((xercesc::BinInputStream* makeStream() const))
{
  const char* uri = std::getenv("OPENMS_TEST_S3_URI");
  if (uri == nullptr || !S3InputSource::isAvailable())
  {
    NOT_TESTABLE
  }
  else
  {
    PeakMap local, remote;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), local);
    MzMLFile().load(uri, remote);
    TEST_EQUAL(remote.size(), local.size())
    TEST_EQUAL(remote.getChromatograms().size(), local.getChromatograms().size())
    TEST_TRUE(remote.getSpectra() == local.getSpectra())
    TEST_EQUAL(S3InputSource::objectSize(uri) > 0, true)
  }
}
END_SECTION

START_SECTION((static Int64 objectSize(const std::string& uri)))
{
  NOT_TESTABLE // tested with makeStream()
}
END_SECTION

START_SECTION((static std::string readRange(const std::string& uri, Int64 offset, Int64 length)))
{
  NOT_TESTABLE // needs a server, see makeStream()
}
END_SECTION

END_TEST
