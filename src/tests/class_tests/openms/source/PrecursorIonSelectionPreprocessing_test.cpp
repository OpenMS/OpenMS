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
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(PrecursorIonSelectionPreprocessing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PrecursorIonSelectionPreprocessing* ptr = nullptr;
PrecursorIonSelectionPreprocessing* nullPointer = nullptr;
START_SECTION(PrecursorIonSelectionPreprocessing())
	ptr = new PrecursorIonSelectionPreprocessing();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PrecursorIonSelectionPreprocessing())
	delete ptr;
END_SECTION

ptr = new PrecursorIonSelectionPreprocessing();

START_SECTION(PrecursorIonSelectionPreprocessing(const PrecursorIonSelectionPreprocessing &source))
	PrecursorIonSelectionPreprocessing copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PrecursorIonSelectionPreprocessing& operator=(const PrecursorIonSelectionPreprocessing &source))
  PrecursorIonSelectionPreprocessing copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

Param param;
param.setValue("precursor_mass_tolerance",0.9);
param.setValue("precursor_mass_tolerance_unit","Da");
param.setValue("missed_cleavages",0);
std::string tmp_filename;
NEW_TMP_FILE(tmp_filename);
param.setValue("preprocessed_db_path",tmp_filename);
ptr->setParameters(param);
ptr->dbPreprocessing(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_db.fasta"),true);

START_SECTION((const std::map<String,std::vector<double> >& getProtMasses() const))
	std::map<String,std::vector<double> > prot_map = ptr->getProtMasses();
	TEST_EQUAL(prot_map.size(), 3)
END_SECTION

START_SECTION((const std::vector<double> & getMasses(String acc) const))
	const std::vector<double>& pep_masses= ptr->getMasses("P01008");
	TEST_EQUAL(pep_masses.size(), 14)
	TEST_REAL_SIMILAR(pep_masses[0],1356.68332791328)
	const std::vector<double>& pep_masses2= ptr->getMasses("P02787");
  TEST_EQUAL(pep_masses2.size(), 19)
	TEST_REAL_SIMILAR(pep_masses2[0],306.159984588623)
END_SECTION


START_SECTION(void dbPreprocessing(String db_path, bool save=true))
	std::map<String,std::vector<double> > prot_map = ptr->getProtMasses();
	TEST_EQUAL(prot_map.size(), 3)
END_SECTION

START_SECTION(double getWeight(double mass))
  double w = ptr->getWeight(147.113);
  TEST_REAL_SIMILAR(w,1)
END_SECTION

START_SECTION(void loadPreprocessing())
	PrecursorIonSelectionPreprocessing ldb;
  ldb.setParameters(param);
  ldb.loadPreprocessing();
  TEST_EQUAL(ldb.getProtMasses().size(),3)
  double w = ldb.getWeight(147.113);
  TEST_REAL_SIMILAR(w,1)

	std::vector<double> pep_masses_l = ldb.getMasses("P01008");
  std::vector<double> pep_masses = ptr->getMasses("P01008");
  TEST_EQUAL(pep_masses_l.size(),pep_masses.size())
	TEST_REAL_SIMILAR(pep_masses_l[0],pep_masses[0])
END_SECTION

PrecursorIonSelectionPreprocessing rt_pt_pp;
rt_pt_pp.setParameters(param);
rt_pt_pp.dbPreprocessing(OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_db.fasta"),
												 OPENMS_GET_TEST_DATA_PATH("PrecursorIonSelectionPreprocessing_rt.model"),
												 OPENMS_GET_TEST_DATA_PATH("DetectabilitySimulation.svm"),false);

START_SECTION(void dbPreprocessing(String db_path,String rt_model_path,String dt_model_path,bool save=true))
  TEST_EQUAL(rt_pt_pp.getProtMasses().size(),3);
  double w = rt_pt_pp.getWeight(147.113);
  TEST_REAL_SIMILAR(w,1)
	TEST_REAL_SIMILAR(rt_pt_pp.getRT("P01008",1),831.46429)
 	TEST_REAL_SIMILAR(rt_pt_pp.getPT("P01008",1),0.0402)
END_SECTION

START_SECTION(double getRT(String prot_id,Size peptide_index))
	TEST_REAL_SIMILAR(rt_pt_pp.getRT("P01008",1),831.46429)
END_SECTION

START_SECTION(double getPT(String prot_id,Size peptide_index))
  TEST_REAL_SIMILAR(rt_pt_pp.getPT("P01008",1),0.0402)
END_SECTION

START_SECTION((const std::map<String, std::vector<double> >& getProteinRTMap() const))
  const std::map<String, std::vector<double> >& rt_map = rt_pt_pp.getProteinRTMap();
  TEST_REAL_SIMILAR(rt_map.find("P01008")->second[1],831.46429)
  TEST_EQUAL(rt_map.size(),3);
END_SECTION

START_SECTION((const std::map<String, std::vector<double> >& getProteinPTMap() const))
  const std::map<String, std::vector<double> >& pt_map = rt_pt_pp.getProteinPTMap();
  TEST_REAL_SIMILAR(pt_map.find("P01008")->second[1],0.0402)
  TEST_EQUAL(pt_map.size(),3);
END_SECTION

START_SECTION((const std::map<String, std::vector<String> >& getProteinPeptideSequenceMap() const))
const std::map<String, std::vector<String> >& map = rt_pt_pp.getProteinPeptideSequenceMap();
   TEST_EQUAL(map.size(),0);
END_SECTION


START_SECTION((void setFixedModifications(StringList & modifications)))
{
  StringList list = ListUtils::create<String>("Carbamidomethylation (C)");
  ptr->setFixedModifications(list);
  const std::map<char, std::vector<String> > & map = ptr->getFixedModifications();
  TEST_EQUAL(map.size(),1)
  TEST_EQUAL(map.begin()->first,'C')
  TEST_EQUAL(map.begin()->second[0],"Carbamidomethylation")
}
END_SECTION

START_SECTION((const std::map<char, std::vector<String> > & getFixedModifications()))
{
  StringList list = ListUtils::create<String>("Oxidation (M)");
  ptr->setFixedModifications(list);
  const std::map<char, std::vector<String> > & map = ptr->getFixedModifications();
  TEST_EQUAL(map.size(),1)
  TEST_EQUAL(map.begin()->first,'M')
  TEST_EQUAL(map.begin()->second[0],"Oxidation")
}
END_SECTION

START_SECTION((void setGaussianParameters(double mu, double sigma)))
{
  ptr->setGaussianParameters(-3.,10.);
  TEST_REAL_SIMILAR(ptr->getGaussMu(),-3.)
  TEST_REAL_SIMILAR(ptr->getGaussSigma(),10.)
}
END_SECTION

START_SECTION((double getGaussMu()))
{
  ptr->setGaussianParameters(-10.,10.);
  TEST_REAL_SIMILAR(ptr->getGaussMu(),-10.)
}
END_SECTION

START_SECTION((double getGaussSigma()))
{
  ptr->setGaussianParameters(-10.,15.);
  TEST_REAL_SIMILAR(ptr->getGaussSigma(),15.)
}
END_SECTION
std::vector< ConvexHull2D > hulls(2);
hulls[0].addPoint(DPosition<2>(810.0,1.0));
hulls[0].addPoint(DPosition<2>(810.0,2.0));
hulls[1].addPoint(DPosition<2>(854.5,1.0));
hulls[1].addPoint(DPosition<2>(854.5,4.0));


START_SECTION((double getRTProbability(String prot_id, Size peptide_index, Feature &feature)))
{
  Feature f;
  f.setRT(831.46);
  f.setConvexHulls(hulls);
  param.setValue("rt_settings:min_rt",800.);
  param.setValue("rt_settings:max_rt",900.);
  param.setValue("rt_settings:rt_step_size",10.);
  rt_pt_pp.setParameters(param);
  rt_pt_pp.setGaussianParameters(0.,1.);
  TEST_REAL_SIMILAR(rt_pt_pp.getRTProbability("P01008",1,f),0.9973)
}
END_SECTION

START_SECTION((double getRTProbability(double pred_rt, Feature &feature)))
{
  Feature f;
  f.setRT(831.46);
  f.setConvexHulls(hulls);
  TEST_REAL_SIMILAR(rt_pt_pp.getRTProbability(831.46429,f),0.9973)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



