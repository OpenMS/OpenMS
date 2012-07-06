// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar$
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProteinResolver, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProteinResolver* ptr = 0;
ProteinResolver* null_ptr = 0;
START_SECTION(ProteinResolver())
{
  ptr = new ProteinResolver();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ProteinResolver())
{
  delete ptr;
}
END_SECTION

START_SECTION((ProteinResolver(const ProteinResolver &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~ProteinResolver()))
{
  // TODO
}
END_SECTION

START_SECTION((ProteinResolver& operator=(const ProteinResolver &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void resolve(std::vector< ISDGroup > &isd_groups, std::vector< MSDGroup > &msd_groups, std::vector< FASTAFile::FASTAEntry > &protein_data, std::vector< Size > &reindexed_proteins, std::vector< Size > &reindexed_peptides, std::vector< ProteinEntry > &protein_nodes, std::vector< PeptideEntry > &peptide_nodes, std::vector< PeptideIdentification > peptide_identifications, ConsensusMap consensus, EnzymaticDigestion digestor, bool id, UInt min_size)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProteinsAndPeptidesmzTab(std::vector< ProteinEntry > &protein_nodes, std::vector< PeptideEntry > &peptide_nodes, std::vector< Size > &reindexed_proteins, std::vector< Size > &reindexed_peptides, std::vector< PeptideIdentification > &peptide_identifications)))
{
  // TODO
}
END_SECTION

START_SECTION((void writePeptideTable(std::vector< PeptideEntry > &peptides, std::vector< Size > &reindexed_peptides, std::vector< PeptideIdentification > &identifications, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void writePeptideTable(std::vector< PeptideEntry > &peptides, std::vector< Size > &reindexed_peptides, ConsensusMap &consensus, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProteinTable(std::vector< ProteinEntry > &proteins, std::vector< Size > &reindexed_proteins, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProteinGroups(std::vector< ISDGroup > &isd_groups, std::vector< MSDGroup > &msd_groups, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void countTargetDecoy(std::vector< MSDGroup > &msd_groups, ConsensusMap &consensus)))
{
  // TODO
}
END_SECTION

START_SECTION((void countTargetDecoy(std::vector< MSDGroup > &msd_groups, std::vector< PeptideIdentification > &peptide_nodes)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
