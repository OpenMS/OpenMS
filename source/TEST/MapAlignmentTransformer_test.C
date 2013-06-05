// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MapAlignmentTransformer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignmentTransformer* ptr = 0;
MapAlignmentTransformer* null_ptr = 0;

TransformationDescription::DataPoints data;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(1.0, 3.0));

TransformationDescription td(data);
Param params;
td.fitModel("linear", params);

const UInt metaIndexRT = MetaInfo::registry().getIndex("RT");

START_SECTION(MapAlignmentTransformer())
{
	ptr = new MapAlignmentTransformer();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MapAlignmentTransformer())
{
	delete ptr;
}
END_SECTION

START_SECTION((static void transformPeakMaps(std::vector< MSExperiment<> > &maps, const std::vector< TransformationDescription > &given_trafos)))
{
  // create experiment
  MSExperiment<> exp;
  MSExperiment<>::SpectrumType spec;

  // first spectrum (MS)
  spec.setRT(11.1);
  spec.setMSLevel(1);
  exp.addSpectrum(spec);

  // second spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(11.5);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  // third spectrum (MS)
  spec.clear(true);
  spec.setRT(12.2);
  spec.setMSLevel(1);
  exp.addSpectrum(spec);

  // forth spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(12.5);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  std::vector<MSExperiment<> > maps;
  maps.push_back(exp);
  maps.push_back(exp);

  std::vector<TransformationDescription> trafos;
  trafos.push_back(td);
  trafos.push_back(td);

  MapAlignmentTransformer::transformPeakMaps(maps, trafos);

  // check the spectra
  TEST_EQUAL(maps[0][0].getRT(), 23.2)
  TEST_EQUAL(maps[0][1].getRT(), 24.0)
  TEST_EQUAL(maps[0][2].getRT(), 25.4)
  TEST_EQUAL(maps[0][3].getRT(), 26.0)

  TEST_EQUAL(maps[1][0].getRT(), 23.2)
  TEST_EQUAL(maps[1][1].getRT(), 24.0)
  TEST_EQUAL(maps[1][2].getRT(), 25.4)
  TEST_EQUAL(maps[1][3].getRT(), 26.0)

  trafos.push_back(td);
  TEST_EXCEPTION(Exception::IllegalArgument, MapAlignmentTransformer::transformPeakMaps(maps, trafos))
}
END_SECTION

START_SECTION((static void transformFeatureMaps(std::vector< FeatureMap<> > &maps, const std::vector< TransformationDescription > &given_trafos)))
{
  FeatureMap<>::FeatureType f;
  FeatureMap<> featMap;

  f.setRT(11.1);
  featMap.push_back(f);

  f.setRT(11.5);
  featMap.push_back(f);

  f.setRT(12.2);
  featMap.push_back(f);

  f.setRT(12.5);
  featMap.push_back(f);

  vector<FeatureMap<> > maps;
  maps.push_back(featMap);
  maps.push_back(featMap);

  std::vector<TransformationDescription> trafos;
  trafos.push_back(td);
  trafos.push_back(td);

  MapAlignmentTransformer::transformFeatureMaps(maps, trafos);

  // check the spectra
  TEST_EQUAL(maps[0][0].getRT(), 23.2)
  TEST_EQUAL(maps[0][1].getRT(), 24.0)
  TEST_EQUAL(maps[0][2].getRT(), 25.4)
  TEST_EQUAL(maps[0][3].getRT(), 26.0)

  TEST_EQUAL(maps[1][0].getRT(), 23.2)
  TEST_EQUAL(maps[1][1].getRT(), 24.0)
  TEST_EQUAL(maps[1][2].getRT(), 25.4)
  TEST_EQUAL(maps[1][3].getRT(), 26.0)

  trafos.push_back(td);
  TEST_EXCEPTION(Exception::IllegalArgument, MapAlignmentTransformer::transformFeatureMaps(maps, trafos))
}
END_SECTION

START_SECTION((static void transformConsensusMaps(std::vector< ConsensusMap > &maps, const std::vector< TransformationDescription > &given_trafos)))
{
  ConsensusFeature cf;
  ConsensusMap consensusMap;

  cf.setRT(11.1);
  consensusMap.push_back(cf);

  cf.setRT(11.5);
  consensusMap.push_back(cf);

  cf.setRT(12.2);
  consensusMap.push_back(cf);

  cf.setRT(12.5);
  consensusMap.push_back(cf);

  vector<ConsensusMap > maps;
  maps.push_back(consensusMap);
  maps.push_back(consensusMap);

  std::vector<TransformationDescription> trafos;
  trafos.push_back(td);
  trafos.push_back(td);

  MapAlignmentTransformer::transformConsensusMaps(maps, trafos);

  // check the spectra
  TEST_EQUAL(maps[0][0].getRT(), 23.2)
  TEST_EQUAL(maps[0][1].getRT(), 24.0)
  TEST_EQUAL(maps[0][2].getRT(), 25.4)
  TEST_EQUAL(maps[0][3].getRT(), 26.0)

  TEST_EQUAL(maps[1][0].getRT(), 23.2)
  TEST_EQUAL(maps[1][1].getRT(), 24.0)
  TEST_EQUAL(maps[1][2].getRT(), 25.4)
  TEST_EQUAL(maps[1][3].getRT(), 26.0)

  trafos.push_back(td);
  TEST_EXCEPTION(Exception::IllegalArgument, MapAlignmentTransformer::transformConsensusMaps(maps, trafos))
}
END_SECTION

START_SECTION((static void transformPeptideIdentifications(std::vector< std::vector< PeptideIdentification > > &maps, const std::vector< TransformationDescription > &given_trafos)))
{
  PeptideIdentification pi;
  std::vector< PeptideIdentification > pIs;

  pi.setMetaValue(metaIndexRT, 11.1);
  pIs.push_back(pi);

  pi.setMetaValue(metaIndexRT, 11.5);
  pIs.push_back(pi);

  pi.setMetaValue(metaIndexRT, 12.2);
  pIs.push_back(pi);

  pi.setMetaValue(metaIndexRT, 12.5);
  pIs.push_back(pi);

  std::vector< std::vector< PeptideIdentification > > maps;
  maps.push_back(pIs);
  maps.push_back(pIs);

  std::vector<TransformationDescription> trafos;
  trafos.push_back(td);
  trafos.push_back(td);

  MapAlignmentTransformer::transformPeptideIdentifications(maps, trafos);

  // check the spectra
  TEST_EQUAL(maps[0][0].getMetaValue(metaIndexRT), 23.2)
  TEST_EQUAL(maps[0][1].getMetaValue(metaIndexRT), 24.0)
  TEST_EQUAL(maps[0][2].getMetaValue(metaIndexRT), 25.4)
  TEST_EQUAL(maps[0][3].getMetaValue(metaIndexRT), 26.0)

  TEST_EQUAL(maps[1][0].getMetaValue(metaIndexRT), 23.2)
  TEST_EQUAL(maps[1][1].getMetaValue(metaIndexRT), 24.0)
  TEST_EQUAL(maps[1][2].getMetaValue(metaIndexRT), 25.4)
  TEST_EQUAL(maps[1][3].getMetaValue(metaIndexRT), 26.0)

  trafos.push_back(td);
  TEST_EXCEPTION(Exception::IllegalArgument, MapAlignmentTransformer::transformPeptideIdentifications(maps, trafos))
}
END_SECTION

START_SECTION((static void transformSinglePeakMap(MSExperiment<> &msexp, const TransformationDescription &trafo)))
{
  MSExperiment<> exp;
  MSExperiment<>::SpectrumType spec;

  // first spectrum (MS)
  spec.setRT(11.1);
  spec.setMSLevel(1);
  exp.addSpectrum(spec);

  // second spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(11.5);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  // third spectrum (MS)
  spec.clear(true);
  spec.setRT(12.2);
  spec.setMSLevel(1);
  exp.addSpectrum(spec);

  // forth spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(12.5);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  MapAlignmentTransformer::transformSinglePeakMap(exp, td);

  // check the spectra
  TEST_EQUAL(exp[0].getRT(), 23.2)
  TEST_EQUAL(exp[1].getRT(), 24.0)
  TEST_EQUAL(exp[2].getRT(), 25.4)
  TEST_EQUAL(exp[3].getRT(), 26.0)

}
END_SECTION

START_SECTION((static void transformSingleFeatureMap(FeatureMap<> &fmap, const TransformationDescription &trafo)))
{
  FeatureMap<>::FeatureType f;
  FeatureMap<> featMap;

  f.setRT(11.1);
  featMap.push_back(f);

  f.setRT(11.5);
  featMap.push_back(f);

  f.setRT(12.2);
  featMap.push_back(f);

  f.setRT(12.5);
  featMap.push_back(f);


  MapAlignmentTransformer::transformSingleFeatureMap(featMap, td);

  // check the spectra
  TEST_EQUAL(featMap[0].getRT(), 23.2)
  TEST_EQUAL(featMap[1].getRT(), 24.0)
  TEST_EQUAL(featMap[2].getRT(), 25.4)
  TEST_EQUAL(featMap[3].getRT(), 26.0)
}
END_SECTION

START_SECTION((static void transformSingleConsensusMap(ConsensusMap &cmap, const TransformationDescription &trafo)))
{
  ConsensusFeature cf;
  ConsensusMap consensusMap;

  cf.setRT(11.1);
  consensusMap.push_back(cf);

  cf.setRT(11.5);
  consensusMap.push_back(cf);

  cf.setRT(12.2);
  consensusMap.push_back(cf);

  cf.setRT(12.5);
  consensusMap.push_back(cf);

  MapAlignmentTransformer::transformSingleConsensusMap(consensusMap, td);

  // check the spectra
  TEST_EQUAL(consensusMap[0].getRT(), 23.2)
  TEST_EQUAL(consensusMap[1].getRT(), 24.0)
  TEST_EQUAL(consensusMap[2].getRT(), 25.4)
  TEST_EQUAL(consensusMap[3].getRT(), 26.0)
}
END_SECTION

START_SECTION((static void transformSinglePeptideIdentification(std::vector< PeptideIdentification > &pepids, const TransformationDescription &trafo)))
{
  PeptideIdentification pi;
  std::vector< PeptideIdentification > pIs;

  pi.setMetaValue(metaIndexRT, 11.1);
  pIs.push_back(pi);

  pi.setMetaValue(metaIndexRT, 11.5);
  pIs.push_back(pi);

  pi.setMetaValue(metaIndexRT, 12.2);
  pIs.push_back(pi);

  pi.setMetaValue(metaIndexRT, 12.5);
  pIs.push_back(pi);

  MapAlignmentTransformer::transformSinglePeptideIdentification(pIs, td);

  // check the spectra
  TEST_EQUAL(pIs[0].getMetaValue(metaIndexRT), 23.2)
  TEST_EQUAL(pIs[1].getMetaValue(metaIndexRT), 24.0)
  TEST_EQUAL(pIs[2].getMetaValue(metaIndexRT), 25.4)
  TEST_EQUAL(pIs[3].getMetaValue(metaIndexRT), 26.0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



