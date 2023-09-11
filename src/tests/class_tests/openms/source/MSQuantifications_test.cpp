// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/MSQuantifications.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSQuantifications, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSQuantifications* ptr = nullptr;
MSQuantifications* null_ptr = nullptr;
START_SECTION(MSQuantifications())
{
	ptr = new MSQuantifications();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MSQuantifications())
{
	delete ptr;
}
END_SECTION

START_SECTION((MSQuantifications(FeatureMap fm, ExperimentalSettings &es, std::vector< DataProcessing > &dps, std::vector< std::vector< std::pair< String, double > > > labels=(std::vector< std::vector< std::pair< String, double > > >()))))
{
  // TODO
}
END_SECTION

START_SECTION((~MSQuantifications()))
{
  // TODO
}
END_SECTION

START_SECTION((MSQuantifications(const MSQuantifications &source)))
{
  // TODO
}
END_SECTION

START_SECTION((MSQuantifications& operator=(const MSQuantifications &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const MSQuantifications &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const MSQuantifications &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<DataProcessing> getDataProcessingList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Assay>& getAssays() const ))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<Assay>& getAssays()))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<ConsensusMap>& getConsensusMaps() const ))
{
  // TODO
}
END_SECTION

START_SECTION((std::vector<ConsensusMap>& getConsensusMaps()))
{
  // TODO
}
END_SECTION

START_SECTION((void setConsensusMaps(const std::vector< ConsensusMap > &)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<FeatureMap >& getFeatureMaps() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const AnalysisSummary& getAnalysisSummary() const ))
{
  // TODO
}
END_SECTION

START_SECTION((AnalysisSummary& getAnalysisSummary()))
{
  // TODO
}
END_SECTION

START_SECTION((void setDataProcessingList(std::vector< DataProcessing > &dpl)))
{
  // TODO
}
END_SECTION

START_SECTION((void setAnalysisSummaryQuantType(QUANT_TYPES r)))
{
  // TODO
}
END_SECTION

START_SECTION((void addConsensusMap(ConsensusMap &m)))
{
  // TODO
}
END_SECTION

START_SECTION((void assignUIDs()))
{
  // TODO
}
END_SECTION

START_SECTION((void registerExperiment(PeakMap &exp, std::vector< std::vector< std::pair< String, double > > > labels)))
{
  // TODO
}
END_SECTION

START_SECTION((void registerExperiment(ExperimentalSettings &es, std::vector< DataProcessing > &dp, std::vector< std::vector< std::pair< String, double > > > labels=(std::vector< std::vector< std::pair< String, double > > >()))))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::AnalysisSummary] AnalysisSummary()))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::AnalysisSummary] AnalysisSummary(const AnalysisSummary &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::AnalysisSummary] virtual ~AnalysisSummary()))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::AnalysisSummary] AnalysisSummary& operator=(const AnalysisSummary &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::Assay] Assay()))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::Assay] Assay(const Assay &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::Assay] virtual ~Assay()))
{
  // TODO
}
END_SECTION

START_SECTION(([MSQuantifications::Assay] Assay& operator=(const Assay &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



