// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ItraqChannelExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqChannelExtractor* ptr = 0;
ItraqChannelExtractor* nullPointer = 0;
START_SECTION(ItraqChannelExtractor())
{
	ptr = new ItraqChannelExtractor();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ItraqChannelExtractor())
{
	delete ptr;
}
END_SECTION

START_SECTION((ItraqChannelExtractor(Int itraq_type)))
{
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX);
	TEST_EQUAL((StringList) ice.getParameters().getValue("channel_active"), StringList::create("114:liver,117:lung"));
  ItraqChannelExtractor ice2(ItraqChannelExtractor::FOURPLEX);
	TEST_EQUAL((StringList) ice2.getParameters().getValue("channel_active"), StringList::create("114:liver,117:lung"));
}
END_SECTION

START_SECTION((ItraqChannelExtractor(Int itraq_type, const Param &param)))
{
	Param p;
	p.setValue("reporter_mass_shift", 0.1234);
	p.setValue("channel_active", StringList::create("121:this is a test"));
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX, p);
	TEST_EQUAL((double) ice.getParameters().getValue("reporter_mass_shift"), 0.1234);
	TEST_EQUAL((StringList) ice.getParameters().getValue("channel_active"), StringList::create("121:this is a test"));
	
	// this should go wrong
	p.setValue("channel_active", StringList::create("120:channel non existent"));	
	TEST_EXCEPTION(Exception::InvalidParameter, ItraqChannelExtractor ice2(ItraqChannelExtractor::EIGHTPLEX, p));	
}
END_SECTION

START_SECTION((ItraqChannelExtractor(const ItraqChannelExtractor &cp)))
{
	Param p;
	p.setValue("reporter_mass_shift", 0.1234);
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX, p);
	ItraqChannelExtractor ice_cp(ice);
	
	TEST_EQUAL(ice_cp.getParameters(), ice.getParameters());
}
END_SECTION

START_SECTION((ItraqChannelExtractor& operator=(const ItraqChannelExtractor &rhs)))
{
	Param p;
	p.setValue("reporter_mass_shift", 0.1234);
  ItraqChannelExtractor ice(ItraqChannelExtractor::EIGHTPLEX, p);
	ItraqChannelExtractor ice_cp;
	ice_cp=ice;
	
	TEST_EQUAL(ice_cp.getParameters(), ice.getParameters());
}
END_SECTION



START_SECTION((void run(const MSExperiment< Peak1D > &ms_exp_data, ConsensusMap &consensus_map)))
{
	MzDataFile mz_data_file;
	MSExperiment<Peak1D > exp;
	mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.mzData"),exp);
	Param p;
	p.setValue("channel_active", StringList::create("114:ref,115:something,116:else"));
	p.setValue("select_activation","");
  ItraqChannelExtractor ice(ItraqChannelExtractor::FOURPLEX, p);
	ConsensusMap cm_out;
	ice.run(exp, cm_out);
	
	ConsensusXMLFile cm_file;
	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);
	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.consensusXML"));
	
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



