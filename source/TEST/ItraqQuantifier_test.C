// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ItraqQuantifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqQuantifier* ptr = 0;
START_SECTION(ItraqQuantifier())
{
	ptr = new ItraqQuantifier();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ItraqQuantifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((ItraqQuantifier(Int itraq_type)))
{
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX);
	TEST_EQUAL((String) iq.getParameters().getValue("isotope_correction")=="true", true);
	TEST_EQUAL((Int) iq.getParameters().getValue("channel_reference"), 113);
  ItraqQuantifier iq2(ItraqQuantifier::FOURPLEX);
	TEST_EQUAL((String) iq2.getParameters().getValue("isotope_correction")=="true", true);
	TEST_EQUAL((Int) iq2.getParameters().getValue("channel_reference"), 114);
}
END_SECTION

START_SECTION((ItraqQuantifier(Int itraq_type, const Param &param)))
{
	Param p;
	p.setValue("isotope_correction_values", StringList::create("114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX, p);
	TEST_EQUAL((StringList) iq.getParameters().getValue("isotope_correction_values"),StringList::create( "114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
	
	// this should go wrong
	p.setValue("isotope_correction_values", StringList::create("114:0/0.3/0 , 116:0.1/0.3/3/0.2"));	
	TEST_EXCEPTION(Exception::InvalidParameter, ItraqQuantifier iq2(ItraqQuantifier::EIGHTPLEX, p));	

}
END_SECTION


START_SECTION((ItraqQuantifier(const ItraqQuantifier &cp)))
{
	Param p;
	p.setValue("isotope_correction_values", StringList::create("114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX, p);

	ItraqQuantifier iq_cp(iq);
	
	TEST_EQUAL(iq_cp.getParameters(), iq.getParameters());

}
END_SECTION

START_SECTION((ItraqQuantifier& operator=(const ItraqQuantifier &rhs)))
{
	Param p;
	p.setValue("isotope_correction_values", StringList::create("114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX, p);

	ItraqQuantifier iq_cp;
	iq_cp = iq;
	
	TEST_EQUAL(iq_cp.getParameters(), iq.getParameters());

}
END_SECTION


START_SECTION((void run(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
  ConsensusXMLFile cm_file;
	ConsensusMap cm_in, cm_out;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.consensusXML"),cm_in);

	std::vector<ProteinIdentification> protein_ids;
	std::vector<PeptideIdentification> peptide_ids;
	String document_id;
	IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.idXML"), protein_ids, peptide_ids, document_id);
	
	ItraqQuantifier iq;
	Param p;
	p.setValue("do_normalization", "true");
	iq.setParameters(p);
	iq.run(cm_in, peptide_ids, protein_ids, cm_out);

	String cm_file_out;// = OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.consensusXML");
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);
	
	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.consensusXML"));
}
END_SECTION

START_SECTION((void run(const ConsensusMap &consensus_map_in, const std::vector< PeptideIdentification > &peptide_ids, const std::vector< ProteinIdentification > &protein_ids, ConsensusMap &consensus_map_out)))
{
	NOT_TESTABLE
	//not implemented yet!
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



