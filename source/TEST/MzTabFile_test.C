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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MzTabFile, "$Id: MzTabFile_test.C 10915 2013-04-04 20:14:57Z aiche $")

/////////////////////////////////////////////////////////////

MzTabFile* ptr = 0;
MzTabFile* null_ptr = 0;
START_SECTION(MzTabFile())
{
	ptr = new MzTabFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(void load(const String& filename, MzTab& mzTab) )
	MzTab mzTab;
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab"), mzTab);
END_SECTION
	
START_SECTION(void store(const String& filename, MzTab& mzTab) )
        {
        // save and reload mzTab
	MzTab mzTab;
	MzTab mzTab_reload;
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab"), mzTab);
	MzTabFile().store(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab_tmp"), mzTab);
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab_tmp"), mzTab_reload);
        }
        {
        // save and reload mzTab
	MzTab mzTab;
	MzTab mzTab_reload;
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_iTRAQ.mzTab"), mzTab);
	MzTabFile().store(OPENMS_GET_TEST_DATA_PATH("MzTabFile_iTRAQ.mzTab_tmp"), mzTab);
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_iTRAQ.mzTab_tmp"), mzTab_reload);
        }
        {
        // save and reload mzTab
	MzTab mzTab;
	MzTab mzTab_reload;
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_merged.mzTab"), mzTab);
	MzTabFile().store(OPENMS_GET_TEST_DATA_PATH("MzTabFile_merged.mzTab_tmp"), mzTab);
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_merged.mzTab_tmp"), mzTab_reload);
	}
	{
	MzTab mzTab;
	MzTab mzTab_reload;
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_opt_columns.mzTab"), mzTab);
	MzTabFile().store(OPENMS_GET_TEST_DATA_PATH("MzTabFile_opt_columns.mzTab_tmp"), mzTab);
	MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_opt_columns.mzTab_tmp"), mzTab_reload);
        }
END_SECTION

START_SECTION(~MzTabFile())
{
	delete ptr;
}
END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



