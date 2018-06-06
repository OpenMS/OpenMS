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
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/ANALYSIS/TARGETED/PSLPFormulation.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PSLPFormulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PSLPFormulation* ptr = nullptr;
PSLPFormulation* nullPointer = nullptr;
START_SECTION(PSLPFormulation())
{
	ptr = new PSLPFormulation();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~PSLPFormulation())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename InputPeakType > void createAndSolveILPForKnownLCMSMapFeatureBased(const FeatureMap &features, const MSExperiment< InputPeakType > &experiment, std::vector< IndexTriple > &variable_indices, std::vector< std::vector< std::pair< Size, Size > > > &mass_ranges, std::set< Int > &charges_set, UInt ms2_spectra_per_rt_bin, std::vector< int > &solution_indices)))
{
	std::set<Int> charges_set;
	charges_set.insert(1);
	
	FeatureMap features;
	PeakMap exp;
	std::vector<PSLPFormulation::IndexTriple > variable_indices;
	std::vector<std::vector<std::pair<Size,Size> > > mass_ranges;
	PSLPFormulation wrapper;
	FeatureMap map;

  std::vector<int> solution_indices;

  // test empty input
	PSLPFormulation wrapper2;
  wrapper2.createAndSolveILPForKnownLCMSMapFeatureBased(features,exp,variable_indices,mass_ranges,charges_set,1,solution_indices);
	TEST_EQUAL(variable_indices.size(),0)
	TEST_EQUAL(solution_indices.size(),0)
	solution_indices.clear();

	// now with the same input as with the offline precursor ion selection (can't test them separately)
	FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_features.featureXML"),map);
	PeakMap raw_data;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OfflinePrecursorIonSelection_raw_data.mzML"),raw_data);
	mass_ranges.clear();
	OfflinePrecursorIonSelection ops;
	ops.getMassRanges(map,raw_data,mass_ranges);
	wrapper.createAndSolveILPForKnownLCMSMapFeatureBased(map,raw_data,variable_indices,mass_ranges,charges_set,1,solution_indices);
	TEST_EQUAL(variable_indices.size(),6)
	variable_indices.clear();
	TEST_EQUAL(solution_indices.size(),3)

}
END_SECTION
		
START_SECTION(([PSLPFormulation::IndexLess] bool operator()(IndexTriple const &left, IndexTriple const &right) const))
{
  PSLPFormulation::IndexTriple a,b;
  a.feature = 1;
  b.feature = 2;

  TEST_EQUAL(PSLPFormulation::IndexLess().operator ()(a,b), true )
  TEST_EQUAL(PSLPFormulation::IndexLess().operator ()(b,a), false )
  TEST_EQUAL(PSLPFormulation::IndexLess().operator ()(a,a), false )
}
END_SECTION

START_SECTION(([PSLPFormulation::ScanLess] bool operator()(IndexTriple const &left, IndexTriple const &right) const))
{
  PSLPFormulation::IndexTriple a,b;
  a.scan = 1;
  b.scan = 2;

  TEST_EQUAL(PSLPFormulation::ScanLess().operator ()(a,b), true )
  TEST_EQUAL(PSLPFormulation::ScanLess().operator ()(b,a), false )
  TEST_EQUAL(PSLPFormulation::ScanLess().operator ()(a,a), false )
}
END_SECTION

START_SECTION(([PSLPFormulation::VariableIndexLess] bool operator()(IndexTriple const &left, IndexTriple const &right) const))
{
  PSLPFormulation::IndexTriple a,b;
  a.variable = 1;
  b.variable = 2;

  TEST_EQUAL(PSLPFormulation::VariableIndexLess().operator ()(a,b), true )
  TEST_EQUAL(PSLPFormulation::VariableIndexLess().operator ()(b,a), false )
  TEST_EQUAL(PSLPFormulation::VariableIndexLess().operator ()(a,a), false )
}
END_SECTION

START_SECTION((void setLPSolver(LPWrapper::SOLVER solver)))
{
  PSLPFormulation lp;
  lp.setLPSolver(LPWrapper::SOLVER_GLPK);
  TEST_EQUAL(lp.getLPSolver(),LPWrapper::SOLVER_GLPK)
}
END_SECTION

START_SECTION((LPWrapper::SOLVER getLPSolver()))
{
  PSLPFormulation lp;
  lp.setLPSolver(LPWrapper::SOLVER_GLPK);
  TEST_EQUAL(lp.getLPSolver(),LPWrapper::SOLVER_GLPK)
}
END_SECTION

START_SECTION((void createAndSolveILPForInclusionListCreation(PrecursorIonSelectionPreprocessing & preprocessing, UInt ms2_spectra_per_rt_bin, UInt max_list_size, FeatureMap & precursors, bool solve_ILP = true)))
{
  Param param;
  param.setValue("precursor_mass_tolerance",0.9);
  param.setValue("precursor_mass_tolerance_unit","Da");
  param.setValue("missed_cleavages",0);
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  param.setValue("preprocessed_db_path",tmp_filename);
  PrecursorIonSelectionPreprocessing rt_pt_pp;
  rt_pt_pp.setParameters(param);
  rt_pt_pp.dbPreprocessing(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_db.fasta"),
                           OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_rt.model"),
                           OPENMS_GET_TEST_DATA_PATH("DetectabilitySimulation.svm"),false);
  FeatureMap precursors;
  PSLPFormulation lp;
  lp.createAndSolveILPForInclusionListCreation(rt_pt_pp, 15, 10, precursors, true);
  TEST_EQUAL(precursors.size(),10)
  TEST_EQUAL(precursors[0].getMetaValue("protein"),"P01008")
  TEST_REAL_SIMILAR(precursors[1].getMZ(),1528.743)
}
END_SECTION

START_SECTION((template < typename InputPeakType > void createAndSolveCombinedLPForKnownLCMSMapFeatureBased(const FeatureMap &features, const MSExperiment< InputPeakType > &experiment, std::vector< IndexTriple > &variable_indices, std::vector< int > &solution_indices, std::vector< std::vector< std::pair< Size, Size > > > &mass_ranges, std::set< Int > &charges_set, UInt ms2_spectra_per_rt_bin, Size step_size=0, bool sequential_order=false)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void updateStepSizeConstraint(Size iteration, UInt step_size)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void updateFeatureILPVariables(FeatureMap & new_features, std::vector<IndexTriple> & variable_indices, std::map<Size,std::vector<String> > & feature_constraints_map)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void updateRTConstraintsForSequentialILP(Size & rt_index, UInt ms2_spectra_per_rt_bin, Size max_rt_index)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void updateCombinedILP(FeatureMap & features, PrecursorIonSelectionPreprocessing & preprocessed_db, std::vector<IndexTriple> & variable_indices, std::vector<String> & new_protein_accs, std::vector<String> & protein_accs, PSProteinInference & prot_inference, Size & variable_counter, std::map<String,std::vector<Size> > & protein_feature_map, Feature& new_feature, std::map<String,Size> & protein_variable_index_map, std::map<String,std::set<String> > & prot_id_counter)))
{
  NOT_TESTABLE
}
END_SECTION
        
START_SECTION((void solveILP(std::vector<int> & solution_indices)))
{
  NOT_TESTABLE
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
