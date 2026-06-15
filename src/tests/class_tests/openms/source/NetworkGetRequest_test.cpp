// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SYSTEM/NetworkGetRequest.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(NetworkGetRequest, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

NetworkGetRequest* ptr = nullptr;
NetworkGetRequest* null_ptr = nullptr;

START_SECTION(NetworkGetRequest())
{
  ptr = new NetworkGetRequest();
  TEST_NOT_EQUAL(ptr, null_ptr)
  // A freshly constructed instance has no error and no response yet.
  TEST_EQUAL(ptr->hasError(), false)
  TEST_EQUAL(ptr->getErrorString().empty(), true)
  TEST_EQUAL(ptr->getResponse().empty(), true)
  TEST_EQUAL(ptr->getResponseBinary().empty(), true)
}
END_SECTION

START_SECTION(~NetworkGetRequest())
{
  delete ptr;
}
END_SECTION

START_SECTION(void setUrl(const std::string& url))
  NOT_TESTABLE // exercised through run() below
END_SECTION

START_SECTION(void setTimeout(int seconds))
  NOT_TESTABLE // exercised through run() below
END_SECTION

START_SECTION(void run())
{
  // Contract (see header): run() never throws. Failures (libcurl init, transport
  // errors, HTTP >= 400) are reported via hasError()/getErrorString(). These cases
  // exercise the offline / bad-URL paths without needing a reachable server.

  // (1) Empty URL: libcurl rejects it (CURLE_URL_MALFORMAT) with no network access.
  {
    NetworkGetRequest req;
    req.setUrl("");
    bool threw = false;
    try { req.run(); } catch (...) { threw = true; }
    TEST_EQUAL(threw, false)
    TEST_EQUAL(req.hasError(), true)
    TEST_NOT_EQUAL(req.getErrorString().size(), 0)
    TEST_EQUAL(req.getResponse().empty(), true)
    TEST_EQUAL(req.getResponseBinary().empty(), true)
  }

  // (2) Malformed URL (scheme but no host): a libcurl URL-parse error, no network.
  {
    NetworkGetRequest req;
    req.setUrl("http://");
    bool threw = false;
    try { req.run(); } catch (...) { threw = true; }
    TEST_EQUAL(threw, false)
    TEST_EQUAL(req.hasError(), true)
    TEST_NOT_EQUAL(req.getErrorString().size(), 0)
  }

  // (3) Unresolvable host: the offline / bad-URL transport path. The ".invalid"
  // TLD is reserved (RFC 6761) to always fail name resolution, so this drives a
  // libcurl transport error (CURLE_COULDNT_RESOLVE_HOST) without depending on any
  // reachable server. The short timeout keeps it fast on slow-DNS hosts.
  {
    NetworkGetRequest req;
    req.setUrl("http://openms-nonexistent-host.invalid/resource");
    req.setTimeout(5);
    bool threw = false;
    try { req.run(); } catch (...) { threw = true; }
    TEST_EQUAL(threw, false)
    TEST_EQUAL(req.hasError(), true)
    TEST_NOT_EQUAL(req.getErrorString().size(), 0)
    TEST_EQUAL(req.getResponse().empty(), true)
  }

  // (4) Reuse: a second run() clears the prior error/response state. After a
  // successful default construction + a failing run, re-running on another bad URL
  // must still report an error (state replaced, not accumulated) and not throw.
  {
    NetworkGetRequest req;
    req.setUrl("");
    req.run();
    TEST_EQUAL(req.hasError(), true)

    req.setUrl("http://another-nonexistent-host.invalid/");
    req.setTimeout(5);
    bool threw = false;
    try { req.run(); } catch (...) { threw = true; }
    TEST_EQUAL(threw, false)
    TEST_EQUAL(req.hasError(), true)
  }
}
END_SECTION

START_SECTION(std::string getResponse() const)
  NOT_TESTABLE // exercised through run() above
END_SECTION

START_SECTION(const std::vector<char>& getResponseBinary() const)
  NOT_TESTABLE // exercised through run() above
END_SECTION

START_SECTION(bool hasError() const)
  NOT_TESTABLE // exercised through run() above
END_SECTION

START_SECTION(std::string getErrorString() const)
  NOT_TESTABLE // exercised through run() above
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
