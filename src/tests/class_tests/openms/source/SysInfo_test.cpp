// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <iostream>

///////////////////////////

using namespace OpenMS;

START_TEST(SysInfo, "$Id$")

START_SECTION(static bool getProcessMemoryConsumption(size_t& mem_virtual))
{
  size_t first, after, final;
  TEST_EQUAL(SysInfo::getProcessMemoryConsumption(first), true);
  std::cout << "Memory consumed initally: " << first << " KB" << std::endl;

  {
    MSExperiment<> exp;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_5_long.mzML"), exp);

    TEST_EQUAL(SysInfo::getProcessMemoryConsumption(after), true);
    std::cout << "Memory consumed after reading 20 MB mzML : " << after << " KB" << std::endl;

    TEST_EQUAL(after - first > 10000, true)
  }

  // just for fun. There is probably no guarantee that we get the whole mem back by the memory manager
  TEST_EQUAL(SysInfo::getProcessMemoryConsumption(final), true);
  std::cout << "Memory consumed after release of MSExperiment: " << final << " KB" << std::endl;

  TEST_EQUAL(after > final, true)

}
END_SECTION

END_TEST
