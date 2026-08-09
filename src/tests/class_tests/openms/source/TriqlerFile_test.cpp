// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lukas Heumos $
// $Authors: Lukas Heumos $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/TriqlerFile.h>

using namespace OpenMS;

START_TEST(TriqlerFile, "$Id$")

START_SECTION(void OpenMS::TriqlerFile::storeLFQ( const std::string &filename,
                                                  const ConsensusMap &consensus_map,
                                                  const OpenMS::ExperimentalDesign& design,
                                                  const StringList& reannotate_filenames,
                                                  const std::string& condition))
{
  // exercised via ProteomicsLFQ -out_triqler (TOPP_ProteomicsLFQ_*), which is the only
  // caller since the TriqlerConverter tool was removed. Note those tests write the Triqler
  // output but never diff it, so no value assertion covers this method yet.
}
END_SECTION

END_TEST
