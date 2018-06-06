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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricQuantifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricQuantifier* ptr = nullptr;
IsobaricQuantifier* null_ptr = nullptr;

ItraqFourPlexQuantitationMethod quant_meth;

START_SECTION((IsobaricQuantifier(const IsobaricQuantitationMethod *const quant_method)))
{
	ptr = new IsobaricQuantifier(&quant_meth);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricQuantifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((IsobaricQuantifier(const IsobaricQuantifier &other)))
{
  IsobaricQuantifier quantifier(&quant_meth);
  IsobaricQuantifier * quantifier2 = new IsobaricQuantifier(quantifier);

	TEST_NOT_EQUAL(quantifier2, null_ptr)
  delete quantifier2;

  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((IsobaricQuantifier& operator=(const IsobaricQuantifier &rhs)))
{
  IsobaricQuantifier quantifier(&quant_meth);
  IsobaricQuantifier quantifier2(&quant_meth);

  quantifier2 = quantifier;

  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void quantify(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
	ConsensusXMLFile cm_file;
	ConsensusMap cm_in, cm_out;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricQuantifier.consensusXML"),cm_in);

	IsobaricQuantifier iq(&quant_meth);
	Param p;
	p.setValue("normalization", "true");
  p.setValue("isotope_correction", "true");
	iq.setParameters(p);
	iq.quantify(cm_in,cm_out);

	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);

	WHITELIST("<?xml-stylesheet");
	// WHITELIST("<?xml-stylesheet,consensusElement id=");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricQuantifier_out.consensusXML"));
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
