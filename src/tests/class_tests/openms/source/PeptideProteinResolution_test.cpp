// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeptideProteinResolution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideProteinResolution* ptr = nullptr;
PeptideProteinResolution* null_ptr = nullptr;
START_SECTION(PeptideProteinResolution())
{
	ptr = new PeptideProteinResolution();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(static void PeptideProteinResolution::run(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides))
{
  vector<ProteinIdentification> prots;
  vector<PeptideIdentification> peps;
  IdXMLFile idf;
  idf.load(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_in.idXML"), prots, peps);  
  PeptideProteinResolution::run(prots, peps);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  IdXMLFile().store(tmp_filename, prots, peps);
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_out.idXML"), tmp_filename);

  prots.clear();
  peps.clear();
  tmp_filename.clear();
  NEW_TMP_FILE(tmp_filename);
  idf.load(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_in2.idXML"), prots, peps);  
  PeptideProteinResolution::run(prots, peps);
  IdXMLFile().store(tmp_filename, prots, peps);
  TEST_FILE_SIMILAR(OPENMS_GET_TEST_DATA_PATH("PeptideProteinResolution_out2.idXML"), tmp_filename);
}
END_SECTION

START_SECTION(~PeptideProteinResolution())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



