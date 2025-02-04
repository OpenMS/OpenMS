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

START_TEST(MSstatsFile, "$Id$")

START_SECTION(void OpenMS::TriqlerFile::storeLFQ( const OpenMS::String &filename, 
                                                  ConsensusMap &consensus_map,
                                                  const OpenMS::ExperimentalDesign& design, 
                                                  const StringList& reannotate_filenames,
                                                  const String& condition,
                                                  const String& retention_time_summarization_method))
{
  // tested via TriqlerConverter tool
}
END_SECTION

END_TEST
