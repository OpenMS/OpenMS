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
// $Maintainer: Hannes Roest, George Rosenberger $
// $Authors: Hannes Roest, George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>

using namespace OpenMS;

//
//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_OpenSwathRTNormalizer OpenSwathRTNormalizer

 @brief The OpenSwathRTNormalizer will find retention time peptides in data.

 This tool will take a description of RT peptides and their normalized
 retention time to write out a transformation file on how to transoform the
 RT space into the normalized space.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_OpenSwathRTNormalizer.cli

 <B>The algorithm parameters for the Analyzer filter are:</B>
 @htmlinclude TOPP_OpenSwathRTNormalizer.html
 
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathRTNormalizer : public TOPPBase
{
public:

  TOPPOpenSwathRTNormalizer() :
  TOPPBase("OpenSwathRTNormalizer",
           "This tool will take a description of RT peptides and their normalized retention time to write out a transformation file on how to transoform the RT space into the normalized space.",
           true)
  {
  }

protected:

  typedef MSExperiment<Peak1D> MapType;

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file with the RT peptides ('TraML' or 'csv')");
    setValidFormats_("tr", ListUtils::create<String>("csv,traML"));
    
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("trafoXML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library)", false);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false);
  }

  void simple_find_best_feature(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
      std::vector<std::pair<double, double> > & pairs)
  {
    for (OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); trgroup_it++)
    {
      // we need at least one feature to find the best one
      OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType * transition_group = &trgroup_it->second;
      if (transition_group->getFeatures().size() == 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Did not find any features for group " + transition_group->getTransitionGroupID());
      }

      MRMFeature * bestf = 0;
      // MRMFeature * secondbest = 0;

      // Find the feature with the highest intensity
      double highest_int = 0;
      for (std::vector<MRMFeature>::iterator mrmfeature = transition_group->getFeaturesMuteable().begin();
           mrmfeature != transition_group->getFeaturesMuteable().end(); mrmfeature++)
      {
        if (mrmfeature->getIntensity() > highest_int)
        {
          bestf = &(*mrmfeature);
          highest_int = mrmfeature->getIntensity();
        }
      }

      // Find the feature with the highest score
      double highest_score = -1000;
      //double second_highest_score = -1000;
      for (std::vector<MRMFeature>::iterator mrmfeature = transition_group->getFeaturesMuteable().begin();
           mrmfeature != transition_group->getFeaturesMuteable().end(); mrmfeature++)
      {
        if (mrmfeature->getOverallQuality() > highest_score)
        {
          //second_highest_score = highest_score;
          bestf = &(*mrmfeature);
          highest_score = mrmfeature->getOverallQuality();
        }
        // else if (mrmfeature->getOverallQuality() > second_highest_score)
        // {
        //   second_highest_score = mrmfeature->getOverallQuality();
        //   secondbest = &(*mrmfeature);
        // }
      }
      String pepref = trgroup_it->second.getTransitions()[0].getPeptideRef();
      if (bestf) {pairs.push_back(std::make_pair(bestf->getRT(), PeptideRTMap[pepref]));}
    }
  }

  ExitCodes main_(int, const char **)
  {

    StringList file_list = getStringList_("in");
    String tr_file_str = getStringOption_("tr");
    String out = getStringOption_("out");
    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");
    const char * tr_file  = tr_file_str.c_str();

    MapType all_xic_maps; // all XICs from all files
    OpenSwath::LightTargetedExperiment targeted_exp;

    std::cout << "Loading TraML file" << std::endl;
    {
      TargetedExperiment transition_exp_;
      TraMLFile().load(tr_file, transition_exp_);
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, targeted_exp);
    }

    // Store the peptide retention times in an intermediate map
    PeptideRTMap.clear();
    for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
    {
      PeptideRTMap[targeted_exp.getPeptides()[i].id] = targeted_exp.getPeptides()[i].rt; 
    }

    MzMLFile f;
    f.setLogType(log_type_);
    TransformationXMLFile trafoxml;
    TransformationDescription trafo;

    // If we have a transformation file, trafo will transform the RT in the
    // scoring according to the model. If we don't have one, it will apply the
    // null transformation.
    if (getStringOption_("rt_norm").size() > 0)
    {
      String trafo_in = getStringOption_("rt_norm");
      String model_type = "linear"; //getStringOption_("model:type");
      trafoxml.load(trafo_in, trafo);
    }

    std::vector<std::pair<double, double> > pairs; // store the RT pairs to write the output trafoXML
    for (Size i = 0; i < file_list.size(); ++i)
    {
      boost::shared_ptr<MapType> swath_map (new MapType()); // the map with the extracted ion chromatograms
      boost::shared_ptr<MapType> xic_map (new MapType());
      FeatureMap<> featureFile;
      std::cout << "RT Normalization working on " << file_list[i] << std::endl;
      f.load(file_list[i], *xic_map.get());

      // get the transitions that we want to use (in swath, only select those
      // from the current window).
      OpenSwath::LightTargetedExperiment transition_exp_used;
      transition_exp_used = targeted_exp;

      std::cout << "nr transitions " << transition_exp_used.getTransitions().size() << std::endl;

      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map;

      MRMFeatureFinderScoring featureFinder;
      Param scoring_params = MRMFeatureFinderScoring().getDefaults();
      scoring_params.setValue("Scores:use_rt_score", "false");
      featureFinder.setParameters(scoring_params);
      
      OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(xic_map);
      featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, trafo, swath_ptr, transition_group_map);

      // add all the chromatograms to the output
      for (Size k = 0; k < xic_map->getChromatograms().size(); k++)
      {
        all_xic_maps.addChromatogram(xic_map->getChromatograms()[k]);
      }

      // find most likely correct feature for each group
      simple_find_best_feature(transition_group_map, pairs);
    }

    std::vector<std::pair<double, double> > pairs_corrected;
    pairs_corrected = MRMRTNormalizer::rm_outliers(pairs, min_rsq, min_coverage, false);

    // store transformation, using a linear model as default
    TransformationDescription trafo_out;
    trafo_out.setDataPoints(pairs_corrected);
    Param model_params;
    model_params.setValue("symmetric_regression", "false");
    String model_type = "linear";
    trafo_out.fitModel(model_type, model_params);
    trafoxml.store(out, trafo_out);

    return EXECUTION_OK;
  }

  std::map<OpenMS::String, double> PeptideRTMap;

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathRTNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
