// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FIAMSScheduler.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FIAMSScheduler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FIAMSScheduler* ptr_1 = nullptr;
FIAMSScheduler* null_ptr_2 = nullptr;
START_SECTION(FIAMSScheduler())
{
    ptr_1 = new FIAMSScheduler(
        String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input/params_test.csv"))
    );
    TEST_NOT_EQUAL(ptr_1, null_ptr_2);
    TEST_EQUAL(ptr_1->getBaseDir(), "/");
}
END_SECTION

START_SECTION(virtual ~FIAMSScheduler())
{
    delete ptr_1;
}
END_SECTION

START_SECTION(FIAMSScheduler)
{
    FIAMSScheduler fia_scheduler(
        String(OPENMS_GET_TEST_DATA_PATH("FIAMS_input/params_test.csv")),
        String(OPENMS_GET_TEST_DATA_PATH(""))
    );
    const vector<map<String, String>> samples = fia_scheduler.getSamples();
    TEST_EQUAL(samples[0].at("time"), "10");
    fia_scheduler.run();
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST