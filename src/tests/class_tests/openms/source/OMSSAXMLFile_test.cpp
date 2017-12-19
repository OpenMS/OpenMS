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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(OMSSAXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

OMSSAXMLFile xml_file;
OMSSAXMLFile* ptr;
OMSSAXMLFile* nullPointer = nullptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications;
vector<PeptideIdentification> peptide_identifications2;
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;

START_SECTION((OMSSAXMLFile()))
	ptr = new OMSSAXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~OMSSAXMLFile())
	delete ptr;
END_SECTION

ptr = new OMSSAXMLFile();

START_SECTION(void setModificationDefinitionsSet(const ModificationDefinitionsSet &rhs))
	ModificationDefinitionsSet mod_set(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)"));
	ptr->setModificationDefinitionsSet(mod_set);
	NOT_TESTABLE
END_SECTION

START_SECTION(void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, bool load_proteins=true, bool load_empty_hits = true))
  // two spectra, first with some hits (mapping to 4 proteins), second is empty
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications);

	TEST_EQUAL(protein_identification.getHits().size(), 4)
	TEST_EQUAL(peptide_identifications.size(), 2)

  xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications, false);
	TEST_EQUAL(protein_identification.getHits().size(), 0)
	TEST_EQUAL(peptide_identifications.size(), 2)

  xml_file.load(OPENMS_GET_TEST_DATA_PATH("OMSSAXMLFile_test_1.xml"),	protein_identification, peptide_identifications, false, false);
	TEST_EQUAL(protein_identification.getHits().size(), 0)
	TEST_EQUAL(peptide_identifications.size(), 1)

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
