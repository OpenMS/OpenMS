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
// $Maintainer: Timo Sachsenberg $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

#include <vector>

///////////////////////////

START_TEST(PTMXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PTMXMLFile* ptr = nullptr;
PTMXMLFile* nullPointer = nullptr;
PTMXMLFile xml_file;

START_SECTION((PTMXMLFile()))
	ptr = new PTMXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((void load(const String& filename, std::map< String, std::pair< String, String > >& ptm_informations)))

	map< String, pair< String, String > > ptm_informations;
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("PTMs.xml"), ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N2O2-CH3")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
END_SECTION


START_SECTION((void store(String filename, std::map< String, std::pair< String, String > > &ptm_informations) const))

	map< String, pair< String, String > > ptm_informations;
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("PTMs.xml"), ptm_informations);
	String temp_filename;
	NEW_TMP_FILE(temp_filename)
	xml_file.store(temp_filename, ptm_informations);
	ptm_informations.clear();
	xml_file.load(temp_filename, ptm_informations);
	
	TEST_EQUAL(ptm_informations["TEST"].first, "N2O2-CH3")
	TEST_EQUAL(ptm_informations["TEST"].second, "KLR")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
