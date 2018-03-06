// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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



