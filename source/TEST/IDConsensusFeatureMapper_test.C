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
#include <OpenMS/ANALYSIS/ID/IDConsensusFeatureMapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
///////////////////////////
	
using namespace OpenMS;
using namespace std;

START_TEST(IDConsensusFeatureMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IDConsensusFeatureMapper* ptr = 0;
CHECK(IDConsensusFeatureMapper())
{
	ptr = new IDConsensusFeatureMapper();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~IDConsensusFeatureMapper())
{
	delete ptr;
}
RESULT

CHECK(void annotate(ConsensusMap& cm,
										const std::vector<PeptideIdentification>& ids,
										const std::vector<ProteinIdentification>& protein_ids,
										CoordinateType mz_delta=0.05,
										CoordinateType rt_delta=0.5,
										bool measure_from_subelements=false))
{
	IDConsensusFeatureMapper mapper;
	FuzzyStringComparator fsc;
	fsc.setAcceptableAbsolute(0.01);
		
	std::vector<ProteinIdentification> protein_ids;
	std::vector<PeptideIdentification> peptide_ids;
	IdXMLFile().load("data/IDConsensusFeatureMapper_in.idXML", protein_ids, peptide_ids);
	
  ConsensusXMLFile cons_file;
  
	{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
	ConsensusMap cons_map;
	cons_file.load("data/IDConsensusFeatureMapper_in.consensusXML", cons_map);
	mapper.annotate(cons_map, peptide_ids, protein_ids);
	cons_file.store(tmp_filename,cons_map);
	TEST_EQUAL(fsc.compare_files(tmp_filename,"data/IDConsensusFeatureMapper_out1.consensusXML"), true);
	}

	{
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
	ConsensusMap cons_map;
	cons_file.load("data/IDConsensusFeatureMapper_in.consensusXML", cons_map);
	mapper.annotate(cons_map, peptide_ids, protein_ids, 0.1, 0.5, true);
	cons_file.store(tmp_filename,cons_map);
	TEST_EQUAL(fsc.compare_files(tmp_filename,"data/IDConsensusFeatureMapper_out2.consensusXML"), true);
	}
	
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



