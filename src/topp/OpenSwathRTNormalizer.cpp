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

    registerFlag_("estimateBestPeptides", "Whether the algorithms should choose the best peptides for normalization. Use this option you do not expect all your peptides to be detected in a sample (e.g. due to them being endogenous peptides or using a less curated list of peptides).", false);

    registerSubsection_("algorithm", "Algorithm parameters section");

    registerSubsection_("peptideEstimation", "Parameters for the peptide estimation (use -estimateBestPeptides to enable).");

    registerSubsection_("outlierDetection", "Parameters for the outlierDetection (use -outlierDetection to enable).");
  }

  Param getSubsectionDefaults_(const String & section) const
  {
    if (section == "algorithm")
    {
      return MRMFeatureFinderScoring().getDefaults();
    }
    else if (section == "peptideEstimation")
    {
      Param p;
      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (section == "outlierDetection")
    {
      Param p;
      p.setValue("useIterativeJackknife", "false", "Whether to use a Jackknife approach optimizing for maximum r-squared when testing for a single outlier (using this option is more expensive with lots of peptides).");
      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when testing for a single outlier (using this option is more stringent and can lead to outliers being retained).");
      p.setValue("useRANSAC", "false", "Whether to use the RANSAC outlier detection algorithm.");
      p.setValue("useRANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      return p;
    }
    return Param();
  }

  void simple_find_best_feature(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
      std::vector<std::pair<double, double> > & pairs, bool useQualCutoff = false, double qualCutoff = 0.0)
  {
    for (OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin(); trgroup_it != transition_group_map.end(); trgroup_it++)
    {
      // we need at least one feature to find the best one
      OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType * transition_group = &trgroup_it->second;
      if (transition_group->getFeatures().size() == 0)
      {
        std::cout << "Did not find any features for group " + transition_group->getTransitionGroupID() << std::endl;
        continue;
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

      if (useQualCutoff && bestf->getOverallQuality() < qualCutoff) {continue;}

      std::cout << " found best feature : " << bestf->getMetaValue("initialPeakQuality") << " and final " << bestf->getOverallQuality() << std::endl;
      String pepref = trgroup_it->second.getTransitions()[0].getPeptideRef();
      if (bestf) {pairs.push_back(std::make_pair(bestf->getRT(), PeptideRTMap[pepref]));}
    }
  }

  bool computeBinnedCoverage(const std::pair<double,double> & rtRange, 
      const std::vector<std::pair<double, double> > & pairs, int nrBins, 
      int minPeptidesPerBin, int minBinsFilled)
  {
    std::vector<int> binCounter(nrBins, 0);
    for (std::vector<std::pair<double, double> >::const_iterator pair_it = pairs.begin(); pair_it != pairs.end(); pair_it++)
    {
      double normRT = (pair_it->second - rtRange.first) / (rtRange.second - rtRange.first); // compute a value between [0,1)
      normRT *= nrBins;
      int bin = (int)normRT;
      if (bin >= nrBins)
      {
        // this should never happen, but just to make sure
        std::cerr << "MRMRTNormalizer::countPeptidesInBins : computed bin was too large (" << bin << "), setting it to the maximum of " << nrBins << std::endl;
        bin = nrBins - 1;
      }
      binCounter[ bin ]++;
    }
    int binsFilled = 0;
    for (Size i = 0; i < binCounter.size(); i++)
    {
      std::cout <<" In bin " << i << " out of " << binCounter.size() << " we have " << binCounter[i] << " peptides " << std::endl;
      if (binCounter[i] >= minPeptidesPerBin) binsFilled++;
    }
    if (binsFilled >= minBinsFilled) return true;
    else return false;
  }

  std::pair<double,double> estimateRTRange(OpenSwath::LightTargetedExperiment & exp)
  {
    if (exp.getPeptides().empty()) 
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "Input list of targets is empty.");
    }
    double max = exp.getPeptides()[0].rt;
    double min = exp.getPeptides()[0].rt;
    for (Size i = 0; i < exp.getPeptides().size(); i++)
    {
      if (exp.getPeptides()[i].rt < min) min = exp.getPeptides()[i].rt;
      if (exp.getPeptides()[i].rt > max) max = exp.getPeptides()[i].rt;
    }
    return std::make_pair(min,max);

  }

  ExitCodes main_(int, const char **)
  {

    StringList file_list = getStringList_("in");
    String tr_file_str = getStringOption_("tr");
    String out = getStringOption_("out");
    DoubleReal min_rsq = getDoubleOption_("min_rsq");
    DoubleReal min_coverage = getDoubleOption_("min_coverage");
    bool estimateBestPeptides = getFlag_("estimateBestPeptides");
    const char * tr_file  = tr_file_str.c_str();

    MapType all_xic_maps; // all XICs from all files
    OpenSwath::LightTargetedExperiment targeted_exp;

    std::cout << "Loading TraML file" << std::endl;
    {
      TargetedExperiment transition_exp_;
      TraMLFile().load(tr_file, transition_exp_);
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, targeted_exp);
    }
    std::pair<double,double> RTRange = estimateRTRange(targeted_exp);
    std::cout << "Rt range from " << RTRange.first << " to " << RTRange.second << std::endl;

    Param pepEstimationParams = getParam_().copy("peptideEstimation:", true);
    Param outlierDetectionParams = getParam_().copy("outlierDetection:", true);

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
      Param scoring_params = getParam_().copy("algorithm:", true);
      scoring_params.setValue("Scores:use_rt_score", "false");
      scoring_params.setValue("Scores:use_elution_model_score", "false");
      if (estimateBestPeptides)
      {
        scoring_params.setValue("TransitionGroupPicker:compute_peak_quality", "true");
        scoring_params.setValue("TransitionGroupPicker:minimal_quality", pepEstimationParams.getValue("InitialQualityCutoff"));
      }
      featureFinder.setParameters(scoring_params);
      featureFinder.setStrictFlag(false);
      
      OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(xic_map);
      featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, trafo, swath_ptr, transition_group_map);

      // add all the chromatograms to the output
      for (Size k = 0; k < xic_map->getChromatograms().size(); k++)
      {
        all_xic_maps.addChromatogram(xic_map->getChromatograms()[k]);
      }

      // find most likely correct feature for each group
      if (estimateBestPeptides)
      {
        simple_find_best_feature(transition_group_map, pairs, true, pepEstimationParams.getValue("OverallQualityCutoff"));
      }
      else
      {
        simple_find_best_feature(transition_group_map, pairs);
      }
    }

    bool enoughPeptides = computeBinnedCoverage(RTRange, pairs, 10, 
        pepEstimationParams.getValue("MinPeptidesPerBin"),
        pepEstimationParams.getValue("MinBinsFilled") );
    // When estimating the best peptides from the data, we need to fulfill the binned coverage
    if (estimateBestPeptides && !enoughPeptides)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "There were not enough bins with the minimal number of peptides");
    }

    std::vector<std::pair<double, double> > pairs_corrected;
    if (outlierDetectionParams.getValue("useRANSAC") == "true")
    {
      double max_rt_threshold = (RTRange.second - RTRange.first) / 30; // estimate of the maximum deviation from RT that is tollerated. Because 120 min gradient can have 4 min elution shift, we use this ratio to find upper RT threshold.
      pairs_corrected = MRMRTNormalizer::rm_outliers_ransac(pairs, min_rsq, min_coverage, outlierDetectionParams.getValue("useRANSACMaxIterations"), max_rt_threshold);
    }
    else
    {
      bool useIterativeJackknife = false;
      if (outlierDetectionParams.getValue("useIterativeJackknife") == "true")
      {
        useIterativeJackknife = true;
      }
      bool useIterativeChauvenet = false;
      if (outlierDetectionParams.getValue("useIterativeChauvenet") == "true")
      {
        useIterativeChauvenet = true;
      }
      pairs_corrected = MRMRTNormalizer::rm_outliers_iterative(pairs, min_rsq, min_coverage, useIterativeChauvenet, useIterativeJackknife);
    }

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
