// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProteinResolver, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProteinResolver* ptr = nullptr;
ProteinResolver* null_ptr = nullptr;
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
