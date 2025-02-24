// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lukas Heumos $
// $Authors: Lukas Heumos $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/MSstatsFile.h>

using namespace OpenMS;

START_TEST(MSstatsFile, "$Id$")

START_SECTION(void OpenMS::MSstatsFile::storeLFQ(const OpenMS::String &filename, ConsensusMap &consensus_map,
                                                 const OpenMS::ExperimentalDesign& design, const StringList& reannotate_filenames,
                                                 const bool is_isotope_label_type, const String& bioreplicate, const String& condition,
                                                 const String& retention_time_summarization_method))
{
  // tested via MSstatsConverter tool
}
END_SECTION

START_SECTION(void OpenMS::MSstatsFile::storeISO(const OpenMS::String &filename, ConsensusMap &consensus_map,
                                                 const OpenMS::ExperimentalDesign& design, const StringList& reannotate_filenames,
                                                 const String& bioreplicate, const String& condition,
                                                 const String& mixture, const String& retention_time_summarization_method))
{
  // tested via MSstatsConverter tool
}
END_SECTION

END_TEST
