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

#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(XTandemXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XTandemXMLFile* ptr;
XTandemXMLFile* nullPointer = nullptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications;

START_SECTION((XTandemXMLFile()))
	ptr = new XTandemXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~XTandemXMLFile())
	delete ptr;
END_SECTION

XTandemXMLFile xml_file;

START_SECTION(void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, ModificationDefinitionsSet& mod_def_set))
{
	ModificationDefinitionsSet mod_set(ListUtils::create<String>(""), ListUtils::create<String>("Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)"));

	xml_file.load(OPENMS_GET_TEST_DATA_PATH("XTandemXMLFile_test.xml"), protein_identification, peptide_identifications, mod_set);
	TEST_EQUAL(peptide_identifications.size(), 303);
	TEST_EQUAL(protein_identification.getHits().size(), 497);
  // should have picked up the default N-terminal modifications:
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 6);
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 0);

  mod_set.setModifications("", "Carbamidomethyl (C),Oxidation (M),Carboxymethyl (C)");
	xml_file.load(OPENMS_GET_TEST_DATA_PATH("XTandemXMLFile_test_2.xml"), protein_identification, peptide_identifications, mod_set);
	TEST_EQUAL(peptide_identifications.size(), 2);
	TEST_EQUAL(protein_identification.getHits().size(), 21);
  // no additional modifications in this case:
  TEST_EQUAL(mod_set.getNumberOfVariableModifications(), 3);
  TEST_EQUAL(mod_set.getNumberOfFixedModifications(), 0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
