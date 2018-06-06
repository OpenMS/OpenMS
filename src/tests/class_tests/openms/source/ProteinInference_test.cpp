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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProteinInference, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProteinInference* ptr = nullptr;
ProteinInference* nullPointer = nullptr;
START_SECTION(ProteinInference())
{
	ptr = new ProteinInference();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ProteinInference())
{
	delete ptr;
}
END_SECTION

START_SECTION((ProteinInference(const ProteinInference &cp)))
{
	NOT_TESTABLE
	// has no members - this is useless
}
END_SECTION

START_SECTION((ProteinInference& operator=(const ProteinInference &rhs)))
{
	NOT_TESTABLE
	// has no members - this is useless
}
END_SECTION

START_SECTION((void infer(ConsensusMap &consensus_map, const UInt reference_map)))
{
	
  ConsensusXMLFile cm_file;
	ConsensusMap cm;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("ProteinInference.consensusXML"),cm);
	
	// delete quantitative info
	for (size_t i=0; i < cm.getProteinIdentifications()[0].getHits().size(); ++i)
	{
		cm.getProteinIdentifications()[0].getHits()[i].clearMetaInfo();
	}
	
	// this should create the quantitation that were in place before deleting them
	ProteinInference inferrer;
	inferrer.infer(cm, 0);
	
	String cm_file_out;// = OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.consensusXML");
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm);
	
	// TOLERANCE_ABSOLUTE(0.01);
	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("ProteinInference.consensusXML"));
	
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



