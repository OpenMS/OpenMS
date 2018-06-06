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
// $Maintainer: Timo Sachsenberg$
// $Authors: Chris Bielow, Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/MSSim.h>
///////////////////////////

#include <OpenMS/KERNEL/RangeUtils.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SIMULATION/SimTypes.h>

#include <algorithm>

// OpenMP support
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

class FindFeature :
  std::unary_function<Feature, bool>
{
public:
  /**
    @brief Constructor

    @param metavalue MetaValue that needs to be present.
    @param reverse if @p reverse is true, operator() returns true if the metavalue does not exist.
  */
  FindFeature(String sequence, Int charge) :
    sequence_(sequence),
    charge_(charge)
  {}

  inline bool operator()(const Feature & f) const
  {
    String f_sequence = "";
    if (f.getPeptideIdentifications().size() > 0)
    {
      if (f.getPeptideIdentifications()[0].getHits().size() > 0)
      {
        f_sequence = f.getPeptideIdentifications()[0].getHits()[0].getSequence().toString();
      }
    }

    return f.getCharge() == charge_ && f_sequence == sequence_;
  }

private:
  String sequence_;
  Int charge_;
};


class FindConsensusFeature :
  std::unary_function<ConsensusFeature, bool>
{
public:
  /**
    @brief Constructor

    @param metavalue MetaValue that needs to be present.
    @param reverse if @p reverse is true, operator() returns true if the metavalue does not exist.
  */
  FindConsensusFeature(String sequence) :
    sequence_(sequence)
  {}

  inline bool operator()(const ConsensusFeature & f) const
  {
    String f_sequence = "";
    if (f.getPeptideIdentifications().size() > 0)
    {
      if (f.getPeptideIdentifications()[0].getHits().size() > 0)
      {
        f_sequence = f.getPeptideIdentifications()[0].getHits()[0].getSequence().toString();
      }
    }

    return f_sequence == sequence_;
  }

private:
  String sequence_;
};


class SumFormulaValue :
  std::unary_function<Feature, bool>
{
public:
  /**
    @brief Constructor

    @param metavalue MetaValue that needs to be present.
    @param reverse if @p reverse is true, operator() returns true if the metavalue does not exist.
  */
  SumFormulaValue(String expected_value) :
    expected_value_(expected_value),
    meta_value_key_("sum_formula")
  {}

  inline bool operator()(const Feature& f) const
  {
    if (f.metaValueExists(meta_value_key_))
    {
      return f.getMetaValue(meta_value_key_) == expected_value_;
    }
    else
    {
      return false;
    }
  }

protected:
  const String expected_value_;
  const String meta_value_key_;
};

START_TEST(MSSim, "$Id$")

// if OpenMP is available features will be generated in different order
// and with this also the PrecursorIonSelection changes leading to`
// differing numbers of MS2 spectra
#ifdef _OPENMP
omp_set_num_threads(1);
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSim* ptr = nullptr;
MSSim* nullPointer = nullptr;
START_SECTION(MSSim())
{
  ptr = new MSSim();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MSSim())
{
  delete ptr;
}
END_SECTION

// we will use this object throughout the test
MSSim mssim;

START_SECTION((void simulate(const SimRandomNumberGenerator &rnd_gen, SimTypes::SampleChannels &peptides)))
{
  SimTypes::MutableSimRandomNumberGeneratorPtr sim_rnd_ptr(new SimTypes::SimRandomNumberGenerator);
  sim_rnd_ptr->initialize(false, false);

  SimTypes::SampleProteins proteins;

  // create some proteins that we want to simulate
  FASTAFile::FASTAEntry protein1;
  protein1.identifier = "1";
  protein1.description = "test-protein-1";
  protein1.sequence = "MTMDKSELVQKAKLAEQAER";

  MetaInfoInterface meta_protein1;
  meta_protein1.setMetaValue("intensity", 1000.0);

  proteins.push_back(SimTypes::SimProtein(protein1, meta_protein1));

  FASTAFile::FASTAEntry protein2;
  protein2.identifier = "2";
  protein2.description = "test-protein-2";
  protein2.sequence = "MTMDKSEVLQKAKIAEQAER";

  MetaInfoInterface meta_protein2;
  meta_protein2.setMetaValue("intensity", 2000.0);

  proteins.push_back(SimTypes::SimProtein(protein2, meta_protein2));

  SimTypes::SampleChannels channels;
  channels.push_back(proteins);

  // TODO: we have to call get parameters first ??
  Param sim_params = mssim.getParameters();
  // define small RT range
  sim_params.setValue("RT:scan_window:min", 210.0);
  sim_params.setValue("RT:scan_window:max", 462.0);

  sim_params.setValue("RawTandemSignal:status", "precursor");

  mssim.setParameters(sim_params);

  mssim.simulate(sim_rnd_ptr, channels);

  // results of simulate are tested individually in the accessors below
  NOT_TESTABLE
}
END_SECTION

START_SECTION((SimTypes::MSSimExperiment const& getExperiment() const ))
{
  // test experiment simulated above

  int nr_ms1 = std::count_if(mssim.getExperiment().begin(),
                             mssim.getExperiment().end(),
                             InMSLevelRange<SimTypes::MSSimExperiment::SpectrumType>(ListUtils::create<Int>("1")));

  int nr_ms2 = std::count_if(mssim.getExperiment().begin(),
                             mssim.getExperiment().end(),
                             InMSLevelRange<SimTypes::MSSimExperiment::SpectrumType>(ListUtils::create<Int>("2")));

#if OPENMS_BOOST_VERSION_MINOR < 56
  TEST_EQUAL(mssim.getExperiment().getNrSpectra(), 230)
  TEST_EQUAL(nr_ms1, 127)
  TEST_EQUAL(nr_ms2, 103)

  TEST_EQUAL(nr_ms2 + nr_ms1, mssim.getExperiment().getNrSpectra())
#else
  TEST_EQUAL(mssim.getExperiment().getNrSpectra(), 234)
  TEST_EQUAL(nr_ms1, 127)
  TEST_EQUAL(nr_ms2, 107)

  TEST_EQUAL(nr_ms2 + nr_ms1, mssim.getExperiment().getNrSpectra())
#endif

  // test empty case when no simulation was performed
  SimTypes::MSSimExperiment empty_experiment;
  MSSim no_sim;
  TEST_EQUAL(no_sim.getExperiment().getSize(), empty_experiment.getSize())
}
END_SECTION

START_SECTION((SimTypes::FeatureMapSim const& getSimulatedFeatures() const ))
{
#if OPENMS_BOOST_VERSION_MINOR < 56
  TEST_EQUAL(mssim.getSimulatedFeatures().size(), 18)
#else
  TEST_EQUAL(mssim.getSimulatedFeatures().size(), 23)
#endif

  // check if all features are contained as expected
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("AKLAEQAER", 3)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("AKLAEQAER", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("AKLAEQAER", 1)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("MTMDK", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("MTMDK", 1)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SELVQKAK", 3)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SELVQKAK", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SELVQKAK", 1)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SEVLQKAK", 3)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SEVLQKAK", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SEVLQKAK", 1)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SEVLQK", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SEVLQK", 1)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SELVQK", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("SELVQK", 1)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("MTMDKSEVLQK", 3)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("MTMDKSEVLQK", 2)) !=  mssim.getSimulatedFeatures().end(), true)
  TEST_EQUAL(find_if(mssim.getSimulatedFeatures().begin(), mssim.getSimulatedFeatures().end(), FindFeature("MTMDKSEVLQK", 1)) !=  mssim.getSimulatedFeatures().end(), true)
}
END_SECTION

START_SECTION((ConsensusMap& getChargeConsensus() ))
{
#if OPENMS_BOOST_VERSION_MINOR < 56
  TEST_EQUAL(mssim.getChargeConsensus().size(), 7)
#else
  TEST_EQUAL(mssim.getChargeConsensus().size(), 9)
#endif

  ConsensusMap::iterator cm_it;

  // AKLAEQAER -> 3 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("AKLAEQAER"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 3)

  // MTMDK -> 2 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("MTMDK"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 2)

  // MTMDKSEVLQK -> 3 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("MTMDKSEVLQK"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 3)

  // SELVQK -> 2 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("SELVQK"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 2)

  // SELVQKAK -> 3 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("SELVQKAK"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 3)

  // SEVLQK -> 3 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("SEVLQK"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 2)

  // SEVLQKAK -> 3 different charge states
  cm_it = find_if(mssim.getChargeConsensus().begin(), mssim.getChargeConsensus().end(), FindConsensusFeature("SEVLQKAK"));
  TEST_EQUAL(cm_it != mssim.getChargeConsensus().end(), true)
  ABORT_IF(cm_it == mssim.getChargeConsensus().end())
  TEST_EQUAL(cm_it->getFeatures().size(), 3)
}
END_SECTION

START_SECTION((ConsensusMap& getLabelingConsensus() ))
{
  // TODO: we need to add another simulation which also labels
}
END_SECTION

START_SECTION((SimTypes::FeatureMapSim const& getContaminants() const ))
{
  TEST_EQUAL(mssim.getContaminants().size(), 37)

  // check expected contaminants are contaiend
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C10H15N1O2S1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C15H24O1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H26O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C18H15O4P1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C22H43N1O1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H32O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H32O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H26O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C24H44N4O4")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C17H34O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C30H58O5S1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C30H58O5S1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C2H6O9Si1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C17H34O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C17H28O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H26O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C2H6O10Si1")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C24H38O4")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C20H38O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H26O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H50O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C24H46O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C17H34O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C17H28O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C17H34O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C16H26O2")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H48O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C20H38O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C20H38O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H48O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H50O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H48O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C20H38O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H48O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H50O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C20H38O7")) != mssim.getContaminants().end(), true)
  TEST_EQUAL(find_if(mssim.getContaminants().begin(), mssim.getContaminants().end(), SumFormulaValue("C26H48O7")) != mssim.getContaminants().end(), true)

}
END_SECTION

START_SECTION((Param getParameters() const ))
{
  Param sim_params = MSSim().getParameters();
  TEST_EQUAL(sim_params.empty(), false)
}
END_SECTION

START_SECTION((SimTypes::MSSimExperiment const& getPeakMap() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void getMS2Identifications(vector<ProteinIdentification>& proteins, vector<PeptideIdentification>& peptides) const))
{
  vector<ProteinIdentification> proteins;
  vector<PeptideIdentification> peptides;

  mssim.getMS2Identifications(proteins, peptides);

  // all 2 proteins should be covered
  TEST_EQUAL(proteins.size(), 1)
  ABORT_IF(proteins.size() != 1)
  TEST_EQUAL(proteins[0].getHits().size(), 2)
  ABORT_IF(proteins[0].getHits().size() != 2)

  // we should have a peptide hit for each ms2 spectrum
  int nr_ms2 = std::count_if(mssim.getExperiment().begin(),
                             mssim.getExperiment().end(),
                             InMSLevelRange<SimTypes::MSSimExperiment::SpectrumType>(ListUtils::create<Int>("2")));
  TEST_EQUAL(peptides.size(), nr_ms2)

  // we assume that there is at least ms2 spectrum that is a mixture of two peptides
  bool is_mixture = false;

  for(vector<PeptideIdentification>::iterator pep_it = peptides.begin();
      pep_it != peptides.end();
      ++pep_it)
  {
    is_mixture |= pep_it->getHits().size() > 1;

    double score = 0.0;
    for(Size i = 0; i < pep_it->getHits().size(); ++i)
    {
      score += pep_it->getHits()[i].getScore();
    }
    // for each PeptideIdentification the sum of scores should be == 1
    TEST_REAL_SIMILAR(score, 1.0)
  }

  // test if there was at least one mix spectrum
  TEST_EQUAL(is_mixture, true)

  // test empty case when no simulation was performed
  MSSim no_sim;
  vector<ProteinIdentification> no_proteins;
  vector<PeptideIdentification> no_peptides;
  no_sim.getMS2Identifications(no_proteins, no_peptides);

  TEST_EQUAL(no_proteins.empty(), true)
  TEST_EQUAL(no_peptides.empty(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
