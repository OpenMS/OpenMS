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
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(OfflinePrecursorIonSelection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OfflinePrecursorIonSelection* ptr = nullptr;
OfflinePrecursorIonSelection* nullPointer = nullptr;
START_SECTION(OfflinePrecursorIonSelection())
	ptr = new OfflinePrecursorIonSelection();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~OfflinePrecursorIonSelection())
	delete ptr;
END_SECTION

ptr = new OfflinePrecursorIonSelection();
std::vector<PeptideIdentification> pep_ids;
std::vector<ProteinIdentification> prot_ids;
//IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_ids.idXML"),prot_ids,pep_ids);

FeatureMap map;
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_features.featureXML"),map);
PeakMap raw_data;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_raw_data.mzML"),raw_data);


START_SECTION((template < typename InputPeakType > void makePrecursorSelectionForKnownLCMSMap(const FeatureMap &features, const MSExperiment< InputPeakType > &experiment, MSExperiment< InputPeakType > &ms2, std::set< Int > &charges_set, bool feature_based)))
{
	PeakMap ms2;
	std::set<Int> charges_set;
	charges_set.insert(1);
	bool feature_based = true;
	Param param;
	param.setValue("ms2_spectra_per_rt_bin", 1);
  param.setValue("min_mz_peak_distance", 3.0);
	ptr->setParameters(param);
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),3)
	TEST_REAL_SIMILAR(ms2[0].getRT(),45)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),55)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),319.19)
	TEST_REAL_SIMILAR(ms2[2].getRT(),65)
	TEST_REAL_SIMILAR(ms2[2].getPrecursors()[0].getMZ(),478.29)
		
	ms2.clear(true);
	feature_based = false;
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),3)
	TEST_REAL_SIMILAR(ms2[0].getRT(),45)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),55)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[2].getRT(),65)
	TEST_REAL_SIMILAR(ms2[2].getPrecursors()[0].getMZ(),336.14)

	ms2.clear(true);
	feature_based = true;
	param.setValue("exclude_overlapping_peaks","true");
	param.setValue("min_mz_peak_distance", 40.);
	ptr->setParameters(param);
	ptr->makePrecursorSelectionForKnownLCMSMap(map,raw_data,ms2,charges_set,feature_based);
	TEST_EQUAL(ms2.size(),2)
	TEST_REAL_SIMILAR(ms2[0].getRT(),45)
	TEST_REAL_SIMILAR(ms2[0].getPrecursors()[0].getMZ(),336.14)
	TEST_REAL_SIMILAR(ms2[1].getRT(),65)
	TEST_REAL_SIMILAR(ms2[1].getPrecursors()[0].getMZ(),478.29)
		
}
END_SECTION	     

START_SECTION((template < typename InputPeakType > void getMassRanges(const FeatureMap &features, const MSExperiment< InputPeakType > &experiment, std::vector< std::vector< std::pair< Size, Size > > > &indices)))
{
	Param param;
	param.setValue("exclude_overlapping_peaks","false");
	ptr->setParameters(param);
	std::vector<std::vector<std::pair<Size,Size> > >  indices;
  FeatureMap map2;
  map2.push_back(map[1]);
  /// test for empty experiment
  PeakMap empty_map;
  TEST_EXCEPTION(Exception::InvalidSize, ptr->getMassRanges(map,empty_map,indices));
  MSSpectrum spec;
  Peak1D p;
  p.setMZ(337.);
  spec.push_back(p);
  p.setMZ(338.);
  spec.push_back(p);
  p.setMZ(339.);
  spec.push_back(p);
  p.setMZ(478.2);
  spec.push_back(p);
  spec.setRT(44.);
  empty_map.addSpectrum(spec);
  spec.setRT(45.);
  empty_map.addSpectrum(spec);
  spec.setRT(46.);
  empty_map.addSpectrum(spec);
  ptr->getMassRanges(map,empty_map,indices);  // led to a memory leak before
  indices.clear();
	ptr->getMassRanges(map,raw_data,indices);
	TEST_EQUAL(indices.size(),3)
	TEST_EQUAL(indices[0][0].first,0)
	TEST_EQUAL(indices[0][0].second,0)
	TEST_EQUAL(indices[0][1].second,0)
	TEST_EQUAL(indices[1][0].first,1)
	TEST_EQUAL(indices[1][0].second,0)
	TEST_EQUAL(indices[1][1].second,0)	
}
END_SECTION

START_SECTION((void setLPSolver(LPWrapper::SOLVER solver)))
{
#if COINOR_SOLVER==1
  ptr->setLPSolver(LPWrapper::SOLVER_COINOR);
  TEST_EQUAL(ptr->getLPSolver(),LPWrapper::SOLVER_COINOR)
#endif
  ptr->setLPSolver(LPWrapper::SOLVER_GLPK);
  TEST_EQUAL(ptr->getLPSolver(),LPWrapper::SOLVER_GLPK)
}
END_SECTION

START_SECTION((LPWrapper::SOLVER getLPSolver()))
{
  // was tested in previous section
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void createProteinSequenceBasedLPInclusionList(String include, String rt_model_file, String pt_model_file, FeatureMap & precursors)))
{
  String include = OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelection_db.fasta");
  String rt_model = OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_rt.model");
  String pt_model	= OPENMS_GET_TEST_DATA_PATH("DetectabilitySimulation.svm");
  FeatureMap precursors;
  ptr->createProteinSequenceBasedLPInclusionList(include,rt_model,pt_model,precursors);
  TEST_EQUAL(precursors.size(),5)
  TEST_EQUAL(precursors[0].getMetaValue("protein"),"P01008")
  TEST_EQUAL(precursors[1].getMetaValue("protein"),"P02787")
  TEST_EQUAL(precursors[4].getMetaValue("protein"),"P10599")
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



