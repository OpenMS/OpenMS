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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/PepXMLFileMascot.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PepXMLFileMascot, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PepXMLFileMascot* ptr = nullptr;
PepXMLFileMascot* nullPointer = nullptr;
PepXMLFileMascot file;
START_SECTION(PepXMLFileMascot())
	ptr = new PepXMLFileMascot();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PepXMLFileMascot())
	delete ptr;
END_SECTION

START_SECTION(void load(const String& filename,  std::map<String, std::vector<AASequence> >& peptides))
	std::map<String, std::vector<AASequence> > peptides;
	std::map<String, std::vector<AASequence> >::iterator it;
	std::vector<AASequence> temp_sequences;
	String filename = OPENMS_GET_TEST_DATA_PATH("PepXMLFileMascot_test.pepXML");
	
	file.load(filename, peptides);
	it = peptides.begin();
	
	TEST_EQUAL(peptides.size(), 3)

	TEST_EQUAL(it->first, "402.42_14_3688.53")
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 4)
	TEST_EQUAL(temp_sequences[0].toUnmodifiedString(), "LAPSAAEDGAFR")
	TEST_EQUAL(temp_sequences[0].toString(), "LAPSAAEDGAFR")
	TEST_EQUAL(temp_sequences[1].toUnmodifiedString(), "GQPLGQLEAHR")			
	TEST_EQUAL(temp_sequences[1].toString(), "GQPLGQLEAHR")
	TEST_EQUAL(temp_sequences[2].toUnmodifiedString(), "SPAPLVPVRLR")			
	TEST_EQUAL(temp_sequences[2].toString(), "SPAPLVPVRLR")
	TEST_EQUAL(temp_sequences[3].toUnmodifiedString(), "QPYLHPSFSK")
	TEST_EQUAL(temp_sequences[3].toString(), "Q(Deamidated)PYLHPSFSK")
	++it;
	TEST_EQUAL(it->first, "404.875_14_333.442")			
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 7)			
	TEST_EQUAL(temp_sequences[0].toUnmodifiedString(), "AQAVAAEIRER")			
	TEST_EQUAL(temp_sequences[0].toString(), "AQAVAAEIRER")
	TEST_EQUAL(temp_sequences[1].toUnmodifiedString(), "AALNAADIVTVR")			
	TEST_EQUAL(temp_sequences[1].toString(), "AALNAADIVTVR")
	TEST_EQUAL(temp_sequences[2].toUnmodifiedString(), "AEPAAELALEAK")			
	TEST_EQUAL(temp_sequences[2].toString(), "AEPAAELALEAK")
	TEST_EQUAL(temp_sequences[3].toUnmodifiedString(), "LANAASPAITQR")			
	TEST_EQUAL(temp_sequences[3].toString(), "LANAASPAITQ(Deamidated)R")
	TEST_EQUAL(temp_sequences[4].toUnmodifiedString(), "SGHSAVLLQDGK")			
	TEST_EQUAL(temp_sequences[4].toString(), "SGHSAVLLQ(Deamidated)DGK")
	TEST_EQUAL(temp_sequences[5].toUnmodifiedString(), "AGAGIANVQAAIR")			
	TEST_EQUAL(temp_sequences[5].toString(), "AGAGIAN(Deamidated)VQ(Deamidated)AAIR")
	TEST_EQUAL(temp_sequences[6].toUnmodifiedString(), "VTAPAARSAALGK")
	TEST_EQUAL(temp_sequences[6].toString(), "VTAPAARSAALGK")
	++it;
	TEST_EQUAL(it->first, "411.766_14_3724.98")			
	temp_sequences = it->second;
	TEST_EQUAL(temp_sequences.size(), 7)
	TEST_EQUAL(temp_sequences[0].toUnmodifiedString(), "LLAWMGRTER")
	TEST_EQUAL(temp_sequences[0].toString(), "LLAWMGRTER")
	TEST_EQUAL(temp_sequences[1].toUnmodifiedString(), "VLALYRAAQAR")
	TEST_EQUAL(temp_sequences[1].toString(), "VLALYRAAQ(Deamidated)AR")			
	TEST_EQUAL(temp_sequences[2].toUnmodifiedString(), "RTLLMSLTGLK")
	TEST_EQUAL(temp_sequences[2].toString(), "RTLLMSLTGLK")		
	TEST_EQUAL(temp_sequences[3].toUnmodifiedString(), "LLGLSRFGLQK")
	TEST_EQUAL(temp_sequences[3].toString(), "LLGLSRFGLQ(Deamidated)K")			
	TEST_EQUAL(temp_sequences[4].toUnmodifiedString(), "MGGIALLDEIGK")
	TEST_EQUAL(temp_sequences[4].toString(), "M(Oxidation)GGIALLDEIGK")			
	TEST_EQUAL(temp_sequences[5].toUnmodifiedString(), "DQMDNALRIR")
	TEST_EQUAL(temp_sequences[5].toString(), "DQMDN(Deamidated)ALRIR")			
	TEST_EQUAL(temp_sequences[6].toUnmodifiedString(), "QTLAGRMVVQK")
	TEST_EQUAL(temp_sequences[6].toString(), "Q(Deamidated)TLAGRMVVQ(Deamidated)K")	
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

