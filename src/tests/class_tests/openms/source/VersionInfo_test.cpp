// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors:  Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/openms_package_version.h>
#include <OpenMS/DATASTRUCTURES/String.h>


/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(VersionInfo, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(static String getTime())
{
  String t = VersionInfo::getTime();
  NOT_TESTABLE
}
END_SECTION

START_SECTION(static String getVersion() )
{
  TEST_STRING_EQUAL(VersionInfo::getVersion(),String(OPENMS_PACKAGE_VERSION).trim());
}
END_SECTION

START_SECTION((static VersionDetails getVersionStruct()))
{
  VersionInfo::VersionDetails detail;
  detail.version_major = 2;
  detail.version_minor = 3;
  detail.version_patch = 0;
  TEST_EQUAL(VersionInfo::getVersionStruct() == detail, true);
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] VersionDetails()))
{
  VersionInfo::VersionDetails detail;
  TEST_EQUAL(detail == VersionInfo::VersionDetails::EMPTY, true)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] VersionDetails(const VersionDetails &other)))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c(detail);
  TEST_EQUAL(c.version_major, detail.version_major)
  TEST_EQUAL(c.version_minor, detail.version_minor)
  TEST_EQUAL(c.version_patch, detail.version_patch)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] bool operator<(const VersionDetails &rhs) const))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_EQUAL(detail < c, false)
  c.version_patch = 3;
  TEST_EQUAL(detail < c, true)
  c.version_patch = 1;
  TEST_EQUAL(detail < c, false)
  c.version_major = 2;
  TEST_EQUAL(detail < c, true)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] bool operator==(const VersionDetails &rhs) const))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_EQUAL(detail == c, true)
  c.version_patch = 3;
  TEST_EQUAL(detail == c, false)
  c.version_patch = 1;
  TEST_EQUAL(detail == c, false)
  c.version_major = 2;
  TEST_EQUAL(detail == c, false)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] bool operator>(const VersionDetails &rhs) const))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_EQUAL(detail > c, false)
  c.version_patch = 3;
  TEST_EQUAL(detail > c, false)
  c.version_patch = 1;
  TEST_EQUAL(detail > c, true)
  c.version_major = 2;
  TEST_EQUAL(detail > c, false)
}
END_SECTION

START_SECTION(([VersionInfo::VersionDetails] static VersionDetails create(const String &version)))
{
  VersionInfo::VersionDetails detail = VersionInfo::VersionDetails::create("1.9.2");
  VersionInfo::VersionDetails c;
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 2;
  TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("1.9");
  c.version_major = 1;
  c.version_minor = 9;
  c.version_patch = 0;
  TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("1.0");
  c.version_major = 1;
  c.version_minor = 0;
  c.version_patch = 0;
  TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("somestring");
  c.version_major = 0;
  c.version_minor = 0;
  c.version_patch = 0;
  TEST_EQUAL(detail == c, true)

  detail = VersionInfo::VersionDetails::create("1.2a.bla");
  c.version_major = 0;
  c.version_minor = 0;
  c.version_patch = 0;
  TEST_EQUAL(detail == c, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
