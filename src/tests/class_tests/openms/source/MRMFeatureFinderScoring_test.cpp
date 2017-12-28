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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
#include <boost/shared_ptr.hpp>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

typedef std::map<String, MRMTransitionGroup< MSChromatogram, OpenSwath::LightTransition > > TransitionGroupMapType;

START_TEST(MRMFeatureFinderScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFeatureFinderScoring* ptr = nullptr;
MRMFeatureFinderScoring* nullPointer = nullptr;

START_SECTION(MRMFeatureFinderScoring())
{
  ptr = new MRMFeatureFinderScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureFinderScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION(void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap< Feature > &output, OpenSwath::LightTargetedExperiment &transition_exp, TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType &transition_group_map))
{
  MRMFeatureFinderScoring ff;
  // Param picker_param = ff.getDefaults();
  // picker_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "legacy"); // old parameters
  // picker_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", 40.0); // old parameters
  // ff.setParameters(picker_param);
      
  MRMFeature feature;
  FeatureMap featureFile;
  TransformationDescription trafo;
  boost::shared_ptr<PeakMap> swath_map (new PeakMap);
  TransitionGroupMapType transition_group_map;
  MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;

  // Load the chromatograms (mzML) and the meta-information (TraML)
  boost::shared_ptr<PeakMap> exp (new PeakMap);
  OpenSwath::LightTargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *exp);
  {
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
  }

  // Pick features in the experiment
#ifdef USE_SP_INTERFACE
  OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  std::vector< OpenSwath::SwathMap > swath_maps(1);
  swath_maps[0].sptr = swath_ptr;
  ff.pickExperiment(chromatogram_ptr, featureFile, transitions, trafo, swath_maps, transition_group_map);
#else
  ff.pickExperiment(exp, featureFile, transitions, trafo, *swath_map, transition_group_map);
#endif

  // Test the number of features found 
  TEST_EQUAL(transition_group_map.size(), 2)

  ///////////////////////////////////////////////////////////////////////////
  //// Scores for the first group
  transition_group = transition_group_map["tr_gr1"];
  TEST_EQUAL(transition_group.size(), 2)
  TEST_EQUAL(transition_group.getFeatures().size(), 1)
  
  // Look closely at the feature we found in the first group
  feature = transition_group.getFeatures()[0];
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR (feature.getRT(), 3119.092);
  TEST_REAL_SIMILAR (feature.getIntensity(), 3615);

  // feature attributes
  TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3089.42993164062);
  TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3154.53002929688);
  TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 3680.16);

  // feature scores
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 0);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_shape"), 0.9981834605);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_rmsd"), 0.108663236);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_corr"), 1);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_elution_model_fit_score"), 0.9854);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_intensity_score"), 0.971);
  TEST_REAL_SIMILAR(feature.getMetaValue("sn_ratio"), 86.0);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_log_sn_score"), 4.45439541136954);

  TOLERANCE_RELATIVE(1.001);
  TEST_REAL_SIMILAR(feature.getMetaValue("rt_score"), 3118.651968);
  TOLERANCE_ABSOLUTE(0.1);

  ///////////////////////////////////////////////////////////////////////////
  //// Scores for the second group
  transition_group = transition_group_map["tr_gr2"];
  TEST_EQUAL(transition_group.size(), 3)
  TEST_EQUAL(transition_group.getFeatures().size(), 2)
  TEST_EQUAL(featureFile.size(), 3)

  // Look closely at the feature we found in the second group
  feature = transition_group.getFeatures()[0];
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(feature.getRT(), 3119.092);
  TEST_REAL_SIMILAR(feature.getIntensity(), 1077.92);

  // feature attributes
  TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3092.85009765625);
  TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3151.10009765625);
  TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 1610.27);

  // feature scores
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 5.70936);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_shape"), 0.7245);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_rmsd"), 0.43566);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_corr"), -0.784);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_elution_model_fit_score"), 0.902);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_intensity_score"), 0.642);
  TEST_REAL_SIMILAR(feature.getMetaValue("sn_ratio"), 30.18);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_log_sn_score"), 3.40718216971789);

  // test legacy parameters
  {
    
    Param picker_param = ff.getDefaults();
    picker_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "legacy"); // old parameters
    picker_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", 40.0); // old parameters
    ff.setParameters(picker_param);

    transition_group_map.clear();
    featureFile.clear();
        
    // Pick features in the experiment
#ifdef USE_SP_INTERFACE
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
    OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    std::vector< OpenSwath::SwathMap > swath_maps(1);
    swath_maps[0].sptr = swath_ptr;
    ff.pickExperiment(chromatogram_ptr, featureFile, transitions, trafo, swath_maps, transition_group_map);
#else
    ff.pickExperiment(exp, featureFile, transitions, trafo, *swath_map, transition_group_map);
#endif

    ///////////////////////////////////////////////////////////////////////////
    //// Scores for the first group
    transition_group = transition_group_map["tr_gr1"];
    TEST_EQUAL(transition_group.size(), 2)
    TEST_EQUAL(transition_group.getFeatures().size(), 1)
    // Look closely at the feature we found in the first group
    feature = transition_group.getFeatures()[0];
    TOLERANCE_ABSOLUTE(0.1);
    TEST_REAL_SIMILAR (feature.getRT(), 3119.092);
    TEST_REAL_SIMILAR (feature.getIntensity(), 3574.23);

    // feature attributes
    TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3096.28);
    TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3147.68);
    TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 3680.16);

    ///////////////////////////////////////////////////////////////////////////
    //// Scores for the second group
    transition_group = transition_group_map["tr_gr2"];
    TEST_EQUAL(transition_group.size(), 3)
    TEST_EQUAL(transition_group.getFeatures().size(), 2)
    TEST_EQUAL(featureFile.size(), 3)
    // Look closely at the feature we found in the second group
    feature = transition_group.getFeatures()[0];
    TOLERANCE_ABSOLUTE(0.1);
    TEST_REAL_SIMILAR(feature.getRT(), 3119.092);
    TEST_REAL_SIMILAR(feature.getIntensity(), 1034.55);

    // feature attributes
    TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3099.7);
    TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3147.68);
    TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 1610.27);
    TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 2.265);
  }

}
END_SECTION

START_SECTION(void pickExperiment(OpenSwath::SpectrumAccessPtr input, FeatureMap< Feature > &output, OpenSwath::LightTargetedExperiment &transition_exp, TransformationDescription trafo, OpenSwath::SpectrumAccessPtr swath_map, TransitionGroupMapType &transition_group_map))
{ 
  MRMFeatureFinderScoring ff;
  Param ff_param = MRMFeatureFinderScoring().getDefaults();
  Param scores_to_use;
  scores_to_use.setValue("use_uis_scores", "true", "Use UIS scores for peptidoform identification ", ListUtils::create<String>("advanced"));
  scores_to_use.setValidStrings("use_uis_scores", ListUtils::create<String>("true,false"));
  ff_param.insert("Scores:", scores_to_use);
  ff_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "legacy"); // old parameters
  ff_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", 40.0); // old parameters
  ff.setParameters(ff_param);

  MRMFeature feature;
  FeatureMap featureFile;
  TransformationDescription trafo;
  boost::shared_ptr<PeakMap> swath_map (new PeakMap);
  TransitionGroupMapType transition_group_map;

  // Load the chromatograms (mzML) and the meta-information (TraML)
  boost::shared_ptr<PeakMap> exp (new PeakMap);
  OpenSwath::LightTargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *exp);
  { 
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_identification_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
  }

  // Pick features in the experiment
#ifdef USE_SP_INTERFACE
  OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  std::vector< OpenSwath::SwathMap > swath_maps(1);
  swath_maps[0].sptr = swath_ptr;
  ff.pickExperiment(chromatogram_ptr, featureFile, transitions, trafo, swath_maps, transition_group_map);
#else
  ff.pickExperiment(exp, featureFile, transitions, trafo, *swath_map, transition_group_map);
#endif

  // Test the number of features found 
  TEST_EQUAL(transition_group_map.size(), 2)

  ///////////////////////////////////////////////////////////////////////////
  //// Scores for the second group
  MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;
  transition_group = transition_group_map["tr_gr2"];
  TEST_EQUAL(transition_group.size(), 3)
  TEST_EQUAL(transition_group.getFeatures().size(), 2)
  TEST_EQUAL(featureFile.size(), 3)

  // Look closely at the feature we found in the second group
  feature = transition_group.getFeatures()[0];
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(feature.getRT(), 3119.092);
  TEST_REAL_SIMILAR(feature.getIntensity(), 1034.55);

  // feature attributes
  TEST_REAL_SIMILAR(feature.getMetaValue("leftWidth" ), 3099.7);
  TEST_REAL_SIMILAR(feature.getMetaValue("rightWidth"), 3147.68);
  TEST_REAL_SIMILAR(feature.getMetaValue("total_xic"), 1610.27);

  // feature scores
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_coelution"), 2.265);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_xcorr_shape"), 0.7245);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_rmsd"), 0.43566);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_library_corr"), -0.784);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_elution_model_fit_score"), 0.902);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_intensity_score"), 0.642);
  TEST_REAL_SIMILAR(feature.getMetaValue("sn_ratio"), 30.18);
  TEST_REAL_SIMILAR(feature.getMetaValue("var_log_sn_score"), 3.40718216971789);

  // feature identification scores
  TEST_EQUAL(feature.getMetaValue("id_target_transition_names"), "tr5;tr2");
  TEST_EQUAL(feature.getMetaValue("id_target_ind_log_intensity"), "5.03352;7.92704");
  TEST_REAL_SIMILAR(feature.getMetaValue("id_target_num_transitions"), 2);
  TEST_EQUAL(feature.getMetaValue("id_target_ind_xcorr_coelution"), "1;1.66667");
  TEST_EQUAL(feature.getMetaValue("id_target_ind_xcorr_shape"), "0.68631;0.690494");
  TEST_EQUAL(feature.getMetaValue("id_target_ind_log_sn_score"), "1.16692;4.45008");
  TEST_EQUAL(feature.getMetaValue("id_target_ind_isotope_correlation"), "");
  TEST_EQUAL(feature.getMetaValue("id_target_ind_isotope_overlap"), "");
  TEST_EQUAL(feature.getMetaValue("id_target_ind_massdev_score"), "");

}
END_SECTION

START_SECTION(void mapExperimentToTransitionList(OpenSwath::SpectrumAccessPtr input, OpenSwath::LightTargetedExperiment &transition_exp, TransitionGroupMapType &transition_group_map, TransformationDescription trafo, double rt_extraction_window))
{
  MRMFeatureFinderScoring ff;
  MRMFeature feature;
  FeatureMap featureFile;
  TransformationDescription trafo;
  PeakMap swath_map;
  TransitionGroupMapType transition_group_map;
  MRMFeatureFinderScoring::MRMTransitionGroupType transition_group;

  // Load the chromatograms (mzML) and the meta-information (TraML)
  boost::shared_ptr<PeakMap> exp (new PeakMap);
  OpenSwath::LightTargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.mzML"), *exp);
  {
    TargetedExperiment transition_exp_;
    TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("OpenSwath_generic_input.TraML"), transition_exp_);
    OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, transitions);
  }

  // Pick features in the experiment
  ff.prepareProteinPeptideMaps_(transitions);
#ifdef USE_SP_INTERFACE
  OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
  ff.mapExperimentToTransitionList(chromatogram_ptr, transitions, transition_group_map, trafo, -1);
#else
  ff.mapExperimentToTransitionList(exp, transitions, transition_group_map, trafo, -1);
#endif

  // Test the number of features found 
  TEST_EQUAL(transition_group_map.size(), 2)

  ///////////////////////////////////////////////////////////////////////////
  //// The first group
  transition_group = transition_group_map["tr_gr1"];
  TEST_EQUAL(transition_group.size(), 2)
  TEST_EQUAL(transition_group.getTransitions().size(), 2)
  TEST_EQUAL(transition_group.getChromatograms().size(), 2)

  TEST_EQUAL(transition_group.hasChromatogram("tr1"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr2"), true)

  TEST_EQUAL(transition_group.getChromatogram("tr2").getNativeID(), "tr2")
  TEST_EQUAL(transition_group.getTransition("tr2").getNativeID(), "tr2")

  ///////////////////////////////////////////////////////////////////////////
  //// The second group
  transition_group = transition_group_map["tr_gr2"];
  TEST_EQUAL(transition_group.size(), 3)
  TEST_EQUAL(transition_group.getTransitions().size(), 3)
  TEST_EQUAL(transition_group.getChromatograms().size(), 3)

  TEST_EQUAL(transition_group.hasChromatogram("tr3"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr4"), true)
  TEST_EQUAL(transition_group.hasChromatogram("tr5"), true)
 
  TEST_EQUAL(transition_group.getChromatogram("tr5").getNativeID(), "tr5")
  TEST_EQUAL(transition_group.getTransition("tr5").getNativeID(), "tr5")

}
END_SECTION

START_SECTION( void scorePeakgroups(MRMTransitionGroupType& transition_group, TransformationDescription & trafo, OpenSwath::SpectrumAccessPtr swath_map, FeatureMap& output) ) 
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(void prepareProteinPeptideMaps_(OpenSwath::LightTargetedExperiment& transition_exp) )
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setStrictFlag(bool f))
  NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

