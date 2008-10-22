// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ItraqQuantifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqQuantifier* ptr = 0;
CHECK(ItraqQuantifier())
{
	ptr = new ItraqQuantifier();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~ItraqQuantifier())
{
	delete ptr;
}
RESULT

CHECK((ItraqQuantifier(Int itraq_type)))
{
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX);
	TEST_EQUAL((String) iq.getParameters().getValue("isotope_correction")=="true", true);
	TEST_EQUAL((Int) iq.getParameters().getValue("channel_reference"), 113);
  ItraqQuantifier iq2(ItraqQuantifier::FOURPLEX);
	TEST_EQUAL((String) iq2.getParameters().getValue("isotope_correction")=="true", true);
	TEST_EQUAL((Int) iq2.getParameters().getValue("channel_reference"), 114);
}
RESULT

CHECK((ItraqQuantifier(Int itraq_type, const Param &param)))
{
	Param p;
	p.setValue("isotope_correction_values", "114:0/0.3/4/0 , 116:0.1/0.3/3/0.2");
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX, p);
	TEST_EQUAL((String) iq.getParameters().getValue("isotope_correction_values"), "114:0/0.3/4/0 , 116:0.1/0.3/3/0.2");
	
	// this should go wrong
	p.setValue("isotope_correction_values", "114:0/0.3/0 , 116:0.1/0.3/3/0.2");	
	TEST_EXCEPTION(Exception::InvalidParameter, ItraqQuantifier iq2(ItraqQuantifier::EIGHTPLEX, p));	

}
RESULT

CHECK((void run(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
  ConsensusXMLFile cm_file;
	ConsensusMap cm_in, cm_out;
	cm_file.load("data/ItraqChannelExtractor.consensusXML",cm_in);
	ItraqQuantifier iq;
	iq.run(cm_in, cm_out);

	String cm_file_out;// = "data/ItraqQuantifier.consensusXML";
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);
	
	FuzzyStringComparator fsc;
	TEST_EQUAL(fsc.compare_files(cm_file_out,"data/ItraqQuantifier.consensusXML"), true);
}
RESULT

CHECK((void run(const ConsensusMap &consensus_map_in, const std::vector< PeptideIdentification > &peptide_ids, const std::vector< ProteinIdentification > &protein_ids, ConsensusMap &consensus_map_out)))
{
	NOT_TESTABLE
	//not implemented yet!
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



