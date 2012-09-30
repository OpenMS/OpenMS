// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

// for setup // TODO re-add this
// #include "OpenSwathTestHelper.h"

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h> 

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h> 
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h>

///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

// TODO : use OpenSwathTestHelper for this
// typedef OpenSWATH_Test::RichPeakChromatogram RichPeakChromatogram;
// typedef OpenSWATH_Test::MRMTransitionGroupType MRMTransitionGroupType;
// typedef OpenSWATH_Test::TransitionType TransitionType;
typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
typedef OpenSwath::LightTransition TransitionType;
typedef MRMTransitionGroup<MSSpectrum, ChromatogramPeak, TransitionType> MRMTransitionGroupType;

namespace OpenSWATH_Test
{
void setup_MRMFeatureFinderScoring(MRMTransitionGroupType& transition_group, std::vector< RichPeakChromatogram >& picked_chroms)
{
  // Load the chromatograms (mzML) and the meta-information (TraML)
  PeakMap exp;
  OpenSwath::LightTargetedExperiment transitions;
  //TargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), exp);

  {
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
    //transitions = transition_exp_;
  }

  // add all the transitions to the peakgroup
  transition_group.setTransitionGroupID("mypeptide");
  transition_group.addTransition(transitions.getTransitions()[0], transitions.getTransitions()[2].getNativeID());
  transition_group.addTransition(transitions.getTransitions()[2], transitions.getTransitions()[0].getNativeID());
  transition_group.addTransition(transitions.getTransitions()[3], transitions.getTransitions()[4].getNativeID());

  // add all the chromatograms to the peakgroup
  {
    Chromatogram chromatogram_old = exp.getChromatograms()[1];
    RichPeakChromatogram chromatogram;
    for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
    {   
      ChromatogramPeak peak;
      peak.setMZ(it->getRT());
      peak.setIntensity(it->getIntensity());
      chromatogram.push_back(peak);
    } 
    //cout << "Chromatogram 0 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
    chromatogram.setMetaValue("product_mz", 618.31);
    chromatogram.setNativeID(chromatogram_old.getNativeID());
    transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
  }
  {
    Chromatogram chromatogram_old = exp.getChromatograms()[0];
    RichPeakChromatogram chromatogram;
    for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
    {   
      ChromatogramPeak peak;
      peak.setMZ(it->getRT());
      peak.setIntensity(it->getIntensity());
      chromatogram.push_back(peak);
    } 
    //cout << "Chromatogram 2 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
    chromatogram.setMetaValue("product_mz", 628.45);
    chromatogram.setNativeID(chromatogram_old.getNativeID());
    transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
  }
  {
    Chromatogram chromatogram_old = exp.getChromatograms()[4]; 
    RichPeakChromatogram chromatogram;
    for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
    {   
      ChromatogramPeak peak;
      peak.setMZ(it->getRT());
      peak.setIntensity(it->getIntensity());
      chromatogram.push_back(peak);
    } 
    //cout << "Chromatogram 3 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
    chromatogram.setMetaValue("product_mz", 651.3);
    chromatogram.setNativeID(chromatogram_old.getNativeID());
    transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
  }

  // do peakpicking, e.g. find a peak at 3120 RT / 170 intensity in all the spectra
  for(Size k = 0; k < transition_group.getChromatograms().size(); k++)
  {
    RichPeakChromatogram picked_chrom;
    ChromatogramPeak peak;
    peak.setMZ(3120);
    peak.setIntensity(170);
    picked_chrom.push_back(peak);

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
    picked_chrom.getFloatDataArrays()[0].push_back(1000.0);
    picked_chrom.getFloatDataArrays()[1].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3140.0);

    picked_chroms.push_back(picked_chrom);
  }

}

}

// typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
// typedef OpenSwath::LightTransition TransitionType;
// typedef MRMTransitionGroup<MSSpectrum, ChromatogramPeak, TransitionType> MRMTransitionGroupType;

//  void setup_MRMFeatureFinderScoring(MRMTransitionGroup<SpectrumT, PeakT, TransitionT> & transition_group, std::vector< ChromatogramT >& picked_chroms)
//  {
//
//    // Load the chromatograms (mzML) and the meta-information (TraML)
//    PeakMap exp;
//    //OpenSwath::LightTargetedExperiment transitions;
//    TargetedExperiment transitions;
//    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), exp);
//
//    {
//      TargetedExperiment transition_exp_;
//      TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
//      //convertTargetedExp(transition_exp_, transitions);
//      transitions = transition_exp_;
//    }
//
//    // add all the transitions to the peakgroup
//    transition_group.setTransitionGroupID("mypeptide");
//    transition_group.addTransition(transitions.getTransitions()[0], transitions.getTransitions()[2].getNativeID());
//    transition_group.addTransition(transitions.getTransitions()[2], transitions.getTransitions()[0].getNativeID());
//    transition_group.addTransition(transitions.getTransitions()[3], transitions.getTransitions()[4].getNativeID());
//
//    // add all the chromatograms to the peakgroup
//    {
//      Chromatogram chromatogram_old = exp.getChromatograms()[1];
//      RichPeakChromatogram chromatogram;
//      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
//      {   
//        ChromatogramPeak peak;
//        peak.setMZ(it->getRT());
//        peak.setIntensity(it->getIntensity());
//        chromatogram.push_back(peak);
//      } 
//      //cout << "Chromatogram 0 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
//      chromatogram.setMetaValue("product_mz", 618.31);
//      chromatogram.setNativeID(chromatogram_old.getNativeID());
//      transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
//    }
//    {
//      Chromatogram chromatogram_old = exp.getChromatograms()[0];
//      RichPeakChromatogram chromatogram;
//      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
//      {   
//        ChromatogramPeak peak;
//        peak.setMZ(it->getRT());
//        peak.setIntensity(it->getIntensity());
//        chromatogram.push_back(peak);
//      } 
//      //cout << "Chromatogram 2 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
//      chromatogram.setMetaValue("product_mz", 628.45);
//      chromatogram.setNativeID(chromatogram_old.getNativeID());
//      transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
//    }
//    {
//      Chromatogram chromatogram_old = exp.getChromatograms()[4]; 
//      RichPeakChromatogram chromatogram;
//      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
//      {   
//        ChromatogramPeak peak;
//        peak.setMZ(it->getRT());
//        peak.setIntensity(it->getIntensity());
//        chromatogram.push_back(peak);
//      } 
//      //cout << "Chromatogram 3 has size " << chromatogram.size() << " and native id " << chromatogram_old.getNativeID() << endl;
//      chromatogram.setMetaValue("product_mz", 651.3);
//      chromatogram.setNativeID(chromatogram_old.getNativeID());
//      transition_group.addChromatogram(chromatogram, chromatogram_old.getNativeID());
//    }
//
//    // do peakpicking, e.g. find a peak at 3120 RT / 170 intensity in all the spectra
//    for(Size k = 0; k < transition_group.getChromatograms().size(); k++)
//    {
//      RichPeakChromatogram picked_chrom;
//      ChromatogramPeak peak;
//      peak.setMZ(3120);
//      peak.setIntensity(170);
//      picked_chrom.push_back(peak);
//
//      picked_chrom.getFloatDataArrays().clear();
//      picked_chrom.getFloatDataArrays().resize(3);
//      picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
//      picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
//      picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
//      picked_chrom.getFloatDataArrays()[0].push_back(1000.0);
//      picked_chrom.getFloatDataArrays()[1].push_back(3100.0);
//      picked_chrom.getFloatDataArrays()[2].push_back(3140.0);
//
//      picked_chroms.push_back(picked_chrom);
//    }
//
//  }

void reorder_transitions(std::vector<TransitionType> & transitions, MRMTransitionGroupType transition_group)
{
  {
    TransitionType t = transition_group.getTransition("tr3");
    t.transition_name = "tr1";
    transitions.push_back(t);
  }
  {
    TransitionType t = transition_group.getTransition("tr1");
    t.transition_name = "tr3";
    transitions.push_back(t);
  }
  {
    TransitionType t = transition_group.getTransition("tr5");
    t.transition_name = "tr3";
    transitions.push_back(t);
  }

}

START_TEST(MRMFeatureScoring, "$Id: MRMFeatureScoring.C 8215 2011-03-29 14:18:26Z aiche $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureScoring* ptr = 0;
MRMFeatureScoring* nullPointer = 0;

START_SECTION(MRMFeatureScoring())
{
  ptr = new MRMFeatureScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION((void MRMFeatureScoring::standardize_data(std::vector<double>& data)))
{
// see seperate test
NOT_TESTABLE
}
END_SECTION

START_SECTION((MRMFeatureScoring::XCorrArrayType MRMFeatureScoring::calcxcorr_new(std::vector<double>& data1, std::vector<double>& data2, int maxdelay, int lag)))
{
// see seperate test
NOT_TESTABLE
}
END_SECTION

START_SECTION((MRMFeatureScoring::XCorrArrayType MRMFeatureScoring::normalizedCalcxcorr(std::vector<double>& data1, std::vector<double>& data2, int maxdelay, int lag)))
{
// see seperate test
NOT_TESTABLE
}
END_SECTION

START_SECTION((MRMFeatureScoring::XCorrArrayType MRMFeatureScoring::calcxcorr(std::vector<double>& data1, std::vector<double>& data2, bool normalize)))
{
// see seperate test
NOT_TESTABLE
}
END_SECTION

///////////////////////////////////////////////////////////////////////////

// From here on we use the setup_MRMFeatureFinderScoring function which we test
// first (giving us "real" data)
START_SECTION((virtual void test_setup()))
{
  // Testing that the setup is correct and the transition group is correctly initialized
  MRMTransitionGroupType transition_group;
  std::vector<RichPeakChromatogram> picked_chroms;

  OpenSWATH_Test::setup_MRMFeatureFinderScoring(transition_group, picked_chroms);

  TEST_EQUAL( transition_group.hasChromatogram("some_unknown_transition"), false)
  TEST_EQUAL(transition_group.hasChromatogram("tr1"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr3"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr5"), true)

  TEST_EQUAL(transition_group.hasTransition("some_unknown_transition"), false)
  TEST_EQUAL(transition_group.hasTransition("tr1"), true)
  TEST_EQUAL(transition_group.hasTransition("tr3"), true)
  TEST_EQUAL(transition_group.hasTransition("tr5"), true)
}
END_SECTION

START_SECTION(
(void initializeXCorrMatrix(MRMFeature& mrmfeature, MRMTransitionGroup< SpectrumType, PeakType >& transition_group, bool normalize)))
{
// see seperate test
NOT_TESTABLE
}
END_SECTION

// testing the individual scores that are produced
// calcXcorrCoelutionScore 
// calcXcorrCoelutionScore_weighted
// calcXcorrShape_score
// calcXcorrShape_score_weighted
// calcLibraryScore
// calcRTScore
// calcElutionFitScore
// calcSNScore
START_SECTION((virtual void test_scores()))
{
  MRMTransitionGroupType transition_group;
  TransformationDescription trafo;
  std::vector< RichPeakChromatogram > picked_chroms;

  OpenSWATH_Test::setup_MRMFeatureFinderScoring(transition_group, picked_chroms);

  // create the corresponding mrm feature
  int chr_idx = 0, peak_idx = 0;
  MRMFeature mrmfeature = MRMTransitionGroupPicker().createMRMFeature(transition_group, picked_chroms, chr_idx, peak_idx);
  TEST_REAL_SIMILAR(mrmfeature.getRT(), 3120.0)

  OpenSwath::IMRMFeature * imrmfeature;
  imrmfeature = new MRMFeatureOpenMS(mrmfeature);

  OpenSwath::ITransitionGroup * itransition_group;
  itransition_group = new TransitionGroupOpenMS<MSSpectrum, ChromatogramPeak, TransitionType>(transition_group);

  //initialize the XCorr Matrix
  MRMFeatureScoring mrmscore;
  mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, true);

  // calculate the normalized library intensity (expected value of the intensities)
  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // (arr1 - mean(arr1) ) / std(arr1)
  // (arr2 - mean(arr2) ) / std(arr2)
  static const double arr_lib[] = {0.5,1,0.5};
  std::vector<double> normalized_library_intensity (arr_lib, arr_lib + sizeof(arr_lib) / sizeof(arr_lib[0]) );
  //mrmscore.standardize_data(normalized_library_intensity);
  double sumx = std::accumulate( normalized_library_intensity.begin(), normalized_library_intensity.end(), 0.0 );
  for(Size m =0; m<normalized_library_intensity.size();m++) { normalized_library_intensity[m] /= sumx;}

  TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionScore(), 2.26491106406735)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionScore_weighted(normalized_library_intensity), 1.375)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrShape_score(), 0.757687954406132)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrShape_score_weighted(normalized_library_intensity), 0.7130856895)

  // numpy
  //data1 = array([1,10000,2000])
  //data2 = array([782.380737304688, 58.3845062255859, 58.3845062255859])
  double library_corr, library_rmsd, d1, d2;
  // We have to reorder the transitions to make the tests work
  std::vector<TransitionType> transitions;
  reorder_transitions(transitions, transition_group);

  mrmscore.calcLibraryScore(imrmfeature, transitions, library_corr, library_rmsd, d1, d2);
  TEST_REAL_SIMILAR(library_corr, -0.654591316)
  TEST_REAL_SIMILAR(library_rmsd, 0.5800337593)

  std::map< OpenMS::String, double > PeptideRTMap; 
#if 0
  // get the id, then get the expected and the experimental retention time
  String native_id = transition_group.getChromatograms()[0].getNativeID();
  TransitionType tr = transition_group.getTransition(native_id);
  const PeptideType * pep = PeptideRefMap[tr.getPeptideRef()];
  double experimental_rt = mrmfeature->getFeature(native_id).getRT();
  double normalized_experimental_rt = trafo.apply(experimental_rt);

  double rt_score = mrmscore.calcRTScore(imrmfeature, itransition_group, trafo, PeptideRTMap);
  TEST_REAL_SIMILAR(rt_score, 3120)
#endif

  /*
   * tested elsewhere!
  EmgScoring emgscore;
  double elution_model_fit_score = emgscore.calcElutionFitScore(mrmfeature, transition_group);
  TEST_REAL_SIMILAR(elution_model_fit_score, 0.924365639)
  */

  // TODO s/n score!
  std::vector<OpenSwath::ISignalToNoisePtr> signal_noise_estimators;
  for (Size k = 0; k < transition_group.getChromatograms().size(); k++)
  {
    OpenSwath::ISignalToNoisePtr snptr(new OpenMS::SignalToNoiseOpenMS<ChromatogramPeak>(transition_group.getChromatograms()[k], 200.0, 30 ));
    signal_noise_estimators.push_back(snptr);
  }

  double sn_score = mrmscore.calcSNScore(imrmfeature, signal_noise_estimators);
  TEST_REAL_SIMILAR(sn_score, 30.180082)
}
END_SECTION

// testing the individual DIA (data independent / SWATH) scores that are produced
// dia_isotope_scores
// dia_massdiff_score
// dia_by_ion_score
// set_dia_parameters
START_SECTION((virtual void test_dia_scores()))
{
  MRMTransitionGroupType transition_group;
  MRMTransitionGroupType transition_group_;
  TransformationDescription trafo;
  std::vector< RichPeakChromatogram > picked_chroms;
  MSExperiment<Peak1D> swath_map;

  OpenSWATH_Test::setup_MRMFeatureFinderScoring(transition_group, picked_chroms);
  transition_group_  = transition_group;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), swath_map);

  int putative_charge_state = 1;
  int by_charge_state = 1;

  // create the corresponding mrm feature
  int chr_idx = 0, peak_idx = 0;
  MRMFeature mrmfeature = MRMTransitionGroupPicker().createMRMFeature(transition_group, picked_chroms, chr_idx, peak_idx);
  TEST_REAL_SIMILAR(mrmfeature.getRT(), 3120.0)

  // find spectrum that is closest to the apex of the peak using binary search
  MSSpectrum<Peak1D> OpenMSspectrum = (*swath_map.RTBegin( mrmfeature.getRT() ));

  OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
  for(MSSpectrum<>::iterator it = OpenMSspectrum.begin(); it != OpenMSspectrum.end(); it++)
  {
    mz_array->data.push_back(it->getMZ());
    intensity_array->data.push_back(it->getIntensity());
  }

  // push back mz first, then intensity.
  // TODO annotate which is which
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  binaryDataArrayPtrs.push_back(mz_array);
  binaryDataArrayPtrs.push_back(intensity_array);

  OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
  sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;
  OpenSwath::SpectrumPtr* spectrum = &sptr;

  MRMFeatureScoring mrmscore;
  OpenSwath::DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity -- because our peptide does not match with the testdata :-)

  // calculate the normalized library intensity (expected value of the intensities)
  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // (arr1 - mean(arr1) ) / std(arr1)
  // (arr2 - mean(arr2) ) / std(arr2)
  static const double arr_lib[] = {0.5,1,0.5};
  std::vector<double> normalized_library_intensity (arr_lib, arr_lib + sizeof(arr_lib) / sizeof(arr_lib[0]) );
  //mrmscore.standardize_data(normalized_library_intensity);
  double sumx = std::accumulate( normalized_library_intensity.begin(), normalized_library_intensity.end(), 0.0 );
  for(Size m =0; m<normalized_library_intensity.size();m++) { normalized_library_intensity[m] /= sumx;}
  
  // Isotope correlation / overlap score: Is this peak part of an
  // isotopic pattern or is it the monoisotopic peak in an isotopic
  // pattern?
  OpenSwath::IMRMFeature * imrmfeature;
  imrmfeature = new MRMFeatureOpenMS(mrmfeature);
  // We have to reorder the transitions to make the tests work
  std::vector<TransitionType> transitions;
  reorder_transitions(transitions, transition_group);
  double isotope_corr = 0, isotope_overlap = 0;
  diascoring.dia_isotope_scores(transitions,
    (*spectrum), imrmfeature, putative_charge_state, isotope_corr, isotope_overlap);

  // Mass deviation score
  double ppm_score = 0, ppm_score_weighted = 0;
  diascoring.dia_massdiff_score(transition_group_.getTransitions(),
    (*spectrum), normalized_library_intensity, ppm_score, ppm_score_weighted);

  // Presence of b/y series score
  double bseries_score = 0, yseries_score = 0;
  String sequence = "SYVAWDR";
  OpenMS::AASequence aas = sequence;
  diascoring.dia_by_ion_score( (*spectrum), aas, by_charge_state, bseries_score, yseries_score);

  TEST_REAL_SIMILAR(isotope_corr, 0.285396985960329* transition_group_.getTransitions().size() )
  TEST_REAL_SIMILAR(isotope_corr, 0.856190957880986)
  TEST_REAL_SIMILAR(isotope_overlap, 0.0599970892071724)

  TEST_REAL_SIMILAR(ppm_score, 1.76388919944981)
  TEST_REAL_SIMILAR(ppm_score_weighted, 0.484116946070573)
  TEST_EQUAL(bseries_score, 0)
  TEST_EQUAL(yseries_score, 1)

  // b/y series score with modifications
  bseries_score = 0, yseries_score = 0;
  aas.setModification(1, "Phospho" ); // modify the Y
  diascoring.dia_by_ion_score( (*spectrum), aas, by_charge_state, bseries_score, yseries_score);
  TEST_EQUAL(bseries_score, 0)
  TEST_EQUAL(yseries_score, 1)
}
END_SECTION

START_SECTION ( XCorrArrayType::iterator xcorrArrayGetMaxPeak(XCorrArrayType array);)
{
  // TODO
}
END_SECTION

START_SECTION (void normalize_sum(double x[], int n);)
{
  // TODO
}
END_SECTION

START_SECTION (void setFitterParam(Param& param))
{
  // TODO
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

