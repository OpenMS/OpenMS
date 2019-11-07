// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Julia Thueringer $
// $Authors: Julia Thueringer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
// calculate pearson distance
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
// create binary tree
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
// align maps and generate output
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

using namespace OpenMS;
using namespace Math;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MapAlignerTree MapAlignerTree

  @brief Corrects retention time distortions between maps, using a tree and identifies features?
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerTreeGuided :
        public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerTreeGuided() :
          TOPPMapAlignerBase("MapAlignerTree", "Tree guided correction of retention time distortions between maps.")
  {
  }

private:
  /// Type to store retention times given for individual peptide sequence
  typedef std::map<String, DoubleList> SeqAndRTList;

  class PeptideIdentificationsPearsonDistance
  {
  public:
    // no const SeqAndRTList because want to access entries
    float operator()(SeqAndRTList& map_first, SeqAndRTList& map_second) const
    {
      // create vectors for both maps containing RTs of identical peptide sequences and
      // get union and intercept amount of peptides
      auto pep1_it = map_first.begin();
      auto pep2_it = map_second.begin();
      vector<double> intercept_rts1;
      vector<double> intercept_rts2;
      float union_size = 0.0;
      while (pep1_it != map_first.end() && pep2_it != map_second.end())
      {
        if (pep1_it->first < pep2_it->first)
        {
          ++pep1_it;
        }
        else if (pep2_it->first < pep1_it->first)
        {
          ++pep2_it;
        }
        else
        {
          // TODO : maybe not one entry for list >1, but min(size1, size2) entries with median?
          double med1 = Math::median(pep1_it->second.begin(), pep1_it->second.end(), true);
          intercept_rts1.push_back(med1);
          double med2 = Math::median(pep2_it->second.begin(), pep2_it->second.end(), true);
          intercept_rts2.push_back(med2);
          ++pep1_it;
          ++pep2_it;
        }
        ++union_size;
      }
      Size intercept_size = intercept_rts1.size();

      // pearsonCorrelationCoefficient(rt_map_i, rt_map_j)
      float pearson_val;
      pearson_val = static_cast<float>(pearsonCorrelationCoefficient(intercept_rts1.begin(), intercept_rts1.end(),
                                                                     intercept_rts2.begin(), intercept_rts2.end()));
      if (pearson_val > 1)
        throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);

      return (1 - (pearson_val * intercept_size / union_size));
    }
  }; // end of PeptideIdentificationsPearsonDifference

  // function declarations
  template <typename MapType>
  void loadInputMaps_(vector<MapType>& maps, StringList& ins, FeatureXMLFile& fxml_file);
  static void getPeptideSequences_(const vector<PeptideIdentification>& peptides, SeqAndRTList& peptide_rts, vector<double>& rts_tmp);
  static void extract_seq_and_rt_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt, vector<double>& maps_ranges);
  static void buildTree_(const vector<FeatureMap>& feature_maps, std::vector<BinaryTreeNode>& tree, vector<double>& maps_ranges);
  void treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree, vector<FeatureMap> feature_maps_transformed,
          vector<double> &maps_ranges, const bool& use_internal_ref, const bool& use_ranges, FeatureMap& map_transformed,
          vector<Size>& trafo_order, const Param& model_params, const String& model_type);
  static void computeTransformationsByOrigInFeature_(vector<FeatureMap>& feature_maps, FeatureMap& map_transformed,
          vector<TransformationDescription>& transformations, const vector<Size>& trafo_order, const Param& model_params,
          const String& model_type);
  void computeTransformedFeatureMaps_(vector<FeatureMap>& feature_maps, const vector<TransformationDescription>& transformations);
  void storeFeatureXMLs_(const vector<FeatureMap>& feature_maps, const StringList& out_files, FeatureXMLFile& fxml_file);
  void storeTransformationDescriptions_(const vector<TransformationDescription>& transformations, StringList& trafos);


  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlags_("featureXML", REF_FLEXIBLE);
    registerStringOption_("use_internal_reference", "string", "false", "Internally the map with the largest rt range is determined as reference for the alignment, otherwise centroid.", false);
    setValidStrings_("use_internal_reference", {"true", "false"});
    registerStringOption_("use_ranges", "string", "true", "Use map with larger 10/90 percentiles range as reference.", false);
    setValidStrings_("use_ranges", {"true", "false"});
    registerSubsection_("align_algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  Param getSubsectionDefaults_(const String&  section) const override {
    if (section == "align_algorithm") {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model") {
      return TOPPMapAlignerBase::getModelDefaults("b_spline");
    }
    return Param(); // this shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes ret = checkParameters_();
    if (ret != EXECUTION_OK) return ret;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in_files = getStringList_("in");
    StringList out_files = getStringList_("out");
    StringList out_trafos = getStringList_("trafo_out");
    bool use_internal_ref = getStringOption_("use_internal_reference").toQString() == "true";
    bool use_ranges = getStringOption_("use_ranges").toQString() == "true";

    Param model_params = getParam_().copy("model:", true);
    String model_type = "b_spline";
    model_params = model_params.copy(model_type+":", true);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    Size in_files_size = in_files.size();
    FeatureXMLFile fxml_file;
    // define here because needed to load and store
    FeatureFileOptions param = fxml_file.getOptions();
    // to save memory don't load convex hulls and subordinates
    param.setLoadSubordinates(false);
    param.setLoadConvexHull(false);
    fxml_file.setOptions(param);

    vector<FeatureMap> feature_maps(in_files_size);
    loadInputMaps_(feature_maps, in_files, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    vector<double> maps_ranges(in_files_size);  // to save ranges for alignment (larger rt_range -> reference)
    std::vector<BinaryTreeNode> tree;    // to construct tree with pearson coefficient
    buildTree_(feature_maps, tree, maps_ranges);

    // print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "  Alignment follows tree: " << ca.newickTree(tree) << endl;

    vector<Size> trafo_order;
    FeatureMap map_transformed;
    // alignment uses copy of feature_maps, returning result within map_transformed (to keep original data)
    treeGuidedAlignment_(tree, feature_maps, maps_ranges, use_internal_ref, use_ranges, map_transformed, trafo_order, model_params, model_type);

    //-------------------------------------------------------------
    // generating output
    //-------------------------------------------------------------
    vector<TransformationDescription> transformations(in_files_size); // for trafo_out
    computeTransformationsByOrigInFeature_(feature_maps, map_transformed, transformations, trafo_order, model_params, model_type);
    computeTransformedFeatureMaps_(feature_maps, transformations);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    // store transformed feature_maps
    storeFeatureXMLs_(feature_maps, out_files, fxml_file);

    // store transformations
    storeTransformationDescriptions_(transformations, out_trafos);

    return EXECUTION_OK;
  }
};

template<typename MapType>
void TOPPMapAlignerTreeGuided::loadInputMaps_(vector<MapType>& maps, StringList& ins, FeatureXMLFile& fxml_file)
{
  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, ins.size(), "loading input files");
  for (Size i = 0; i < ins.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.load(ins[i], maps[i]);
  }
  progresslogger.endProgress();
}

void TOPPMapAlignerTreeGuided::getPeptideSequences_(const vector<PeptideIdentification>& peptides,
                                              TOPPMapAlignerTreeGuided::SeqAndRTList& peptide_rts, vector<double>& rts_tmp)
{
  for (const auto & peptide : peptides)
  {
    if (!peptide.getHits().empty())
    {
      const String& sequence = peptide.getHits()[0].getSequence().toString();
      double rt = peptide.getRT();
      peptide_rts[sequence].push_back(rt);
      rts_tmp.push_back(rt);
    }
  }
}

void TOPPMapAlignerTreeGuided::extract_seq_and_rt_(const vector<FeatureMap>& feature_maps,
                                             vector<SeqAndRTList>& maps_seqAndRt, vector<double>& maps_ranges)
{
  for (Size i = 0; i < feature_maps.size(); ++i)
  {
    double percentile10;
    double percentile90;
    vector<double> rts_tmp(feature_maps[i].size());
    for (auto feature_it = feature_maps[i].begin(); feature_maps[i].end() != feature_it; ++feature_it)
    {
      if (!feature_it->getPeptideIdentifications().empty()) {
        getPeptideSequences_(feature_it->getPeptideIdentifications(), maps_seqAndRt[i], rts_tmp);
      }
    }
    sort(rts_tmp.begin(), rts_tmp.end());

    percentile10 = rts_tmp[rts_tmp.size()*0.1];
    percentile90 = rts_tmp[rts_tmp.size()*0.9];

    maps_ranges[i] = percentile90 - percentile10;
    rts_tmp.clear();
  }
}

void TOPPMapAlignerTreeGuided::buildTree_(const vector<FeatureMap>& feature_maps, std::vector<BinaryTreeNode>& tree, vector<double>& maps_ranges)
{
  vector<SeqAndRTList> maps_seqAndRt(feature_maps.size());
  extract_seq_and_rt_(feature_maps, maps_seqAndRt, maps_ranges);
  PeptideIdentificationsPearsonDistance pepDist;
  SingleLinkage sl;
  DistanceMatrix<float> dist_matrix;
  ClusterHierarchical ch;
  ch.cluster<SeqAndRTList, PeptideIdentificationsPearsonDistance>(maps_seqAndRt, pepDist, sl, tree, dist_matrix);
}

void TOPPMapAlignerTreeGuided::treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree, vector<FeatureMap> feature_maps_transformed,
        vector<double> &maps_ranges, const bool& use_internal_ref, const bool& use_ranges, FeatureMap& map_transformed,
        vector<Size>& trafo_order, const Param& model_params, const String& model_type)
{
  vector<TransformationDescription> trafo_tmp; // use to align
  Size last_trafo = 0;  // look up transformation order in map_sets

  vector<TransformationDescription> transformations_align;  // temporary for aligner output
  Size ref;
  Size to_transform;
  // helper to reorganize rt transformations
  vector<vector<Size>> map_sets(feature_maps_transformed.size());
  for (Size i = 0; i < feature_maps_transformed.size(); ++i)
  {
    vector<Size> tmp;
    tmp.push_back(i);
    map_sets[i] = tmp;
  }

  MapAlignmentAlgorithmIdentification algorithm;
  Param algo_params = getParam_().copy("align_algorithm:", true);
  algorithm.setParameters(algo_params);
  algorithm.setLogType(log_type_);

  for (const auto& node : tree)
  {
    // ----------------
    // prepare alignment
    // ----------------
    vector<FeatureMap> to_align;
    if (use_ranges)
    {
      OPENMS_LOG_INFO << "use larger range"  << endl;
      //  determine the map with larger RT range for 10/90 percentile (->reference)
      if (maps_ranges[node.left_child] > maps_ranges[node.right_child]) {
        ref = node.left_child;
        to_transform = node.right_child;
        maps_ranges[node.right_child] = maps_ranges[node.left_child]; // after transformation same range for both maps :(
        // TODO: check, whether new calculation of range is better
      } else {
        ref = node.right_child;
        to_transform = node.left_child;
        maps_ranges[node.left_child] = maps_ranges[node.right_child]; // after transformation same range for both maps
      }
    }
    else{
      OPENMS_LOG_INFO << "use larger map size"  << endl;
      //  determine the map with more features (->reference)
      if (feature_maps_transformed[node.left_child].size() > feature_maps_transformed[node.right_child].size()) {
        ref = node.left_child;
        to_transform = node.right_child;
      } else {
        ref = node.right_child;
        to_transform = node.left_child;
      }
    }

    last_trafo = to_transform;
    to_align.push_back(feature_maps_transformed[to_transform]);
    to_align.push_back(feature_maps_transformed[ref]);

    // ----------------
    // perform alignment
    // ----------------
    if (use_internal_ref)
    {
      OPENMS_LOG_INFO << "use internal reference for alignment"  << endl;
      // with set reference
      algorithm.align(to_align, transformations_align, 1);
      // transform retention times of non-identity for next iteration
      transformations_align[0].fitModel(model_type, model_params);
      // needed for following iteration steps; store_original_rt is default false
      MapAlignmentTransformer::transformRetentionTimes(feature_maps_transformed[to_transform],
                                                       transformations_align[0], true);
    }
    else{
      OPENMS_LOG_INFO << "no use of an internal reference for alignment"  << endl;
      // without set reference
      algorithm.align(to_align, transformations_align);

      // transform retention times of non-identity for next iteration
      transformations_align[0].fitModel(model_type, model_params);
      transformations_align[1].fitModel(model_type, model_params);

      // needed for following iteration steps; store_original_rt is default false
      MapAlignmentTransformer::transformRetentionTimes(feature_maps_transformed[to_transform],
                                                       transformations_align[0], true);
      MapAlignmentTransformer::transformRetentionTimes(feature_maps_transformed[ref],
                                                       transformations_align[1], true);
    }

    // combine aligned maps, store in both, because tree always calls smaller number
    // also possible: feature_maps_transformed[smallerNumber] = ..[ref]+..[to_transform]
    // or use pointer
    feature_maps_transformed[ref] += feature_maps_transformed[to_transform];
    feature_maps_transformed[ref].updateRanges();
    feature_maps_transformed[to_transform] = feature_maps_transformed[ref];

    // update transformation order for each map
    map_sets[ref].insert(map_sets[ref].end(), map_sets[to_transform].begin(), map_sets[to_transform].end());
    map_sets[to_transform] = map_sets[ref];

    transformations_align.clear();
    to_align.clear();
  }
  // copy last transformed FeatureMap for reference return
  map_transformed = feature_maps_transformed[last_trafo];
  trafo_order = map_sets[last_trafo];
}

void TOPPMapAlignerTreeGuided::computeTransformationsByOrigInFeature_(vector<FeatureMap>& feature_maps, FeatureMap& map_transformed,
                                                                     vector<TransformationDescription>& transformations, const vector<Size>& trafo_order, const Param& model_params,
                                                                     const String& model_type)
{
  FeatureMap::const_iterator fit = map_transformed.begin();
  for (auto & map_idx : trafo_order)
  {
    TransformationDescription::DataPoints trafo_data_tmp;

    for (Size i = 0; i < feature_maps[map_idx].size(); ++i)
    {
      TransformationDescription::DataPoint point;
      if (fit->metaValueExists("original_RT"))
      {
        point.first = fit->getMetaValue("original_RT");
      }
      else{
        point.first = fit->getRT();
      }
      point.second = fit->getRT();
      point.note = fit->getUniqueId();
      trafo_data_tmp.push_back(point);
      ++fit;
    }
    transformations[map_idx] = TransformationDescription(trafo_data_tmp);
    transformations[map_idx].fitModel(model_type, model_params);
    trafo_data_tmp.clear();
  }
}

void TOPPMapAlignerTreeGuided::computeTransformedFeatureMaps_(vector<FeatureMap>& feature_maps, const vector<TransformationDescription>& transformations) {

  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, feature_maps.size(), "computing feature maps");
  for (Size i = 0; i < feature_maps.size(); ++i)
  {
    progresslogger.setProgress(i);
    MapAlignmentTransformer::transformRetentionTimes(feature_maps[i], transformations[i], true);
    //feature_maps[i].updateRanges(); // without: LeakSanitizer detects memory leaks
    // annotate output with data processing info
    addDataProcessing_(feature_maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));
  }
  progresslogger.endProgress();
}

void TOPPMapAlignerTreeGuided::storeFeatureXMLs_(const vector<FeatureMap>& feature_maps, const StringList& out_files, FeatureXMLFile& fxml_file) {

  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, feature_maps.size(), "writing output files");
  for (Size i = 0; i < out_files.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.store(out_files[i], feature_maps[i]);
  }
  progresslogger.endProgress();
}

void TOPPMapAlignerTreeGuided::storeTransformationDescriptions_(const vector<TransformationDescription>& transformations,
                                                          StringList& trafos) {
  // custom progress logger for this task:
  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
  progresslogger.startProgress(0, trafos.size(),
                               "writing transformation files");
  for (Size i = 0; i < trafos.size(); ++i)
  {
    progresslogger.setProgress(i);
    TransformationXMLFile().store(trafos[i], transformations[i]);
  }
  progresslogger.endProgress();
}


int main(int argc, const char** argv)
{
  TOPPMapAlignerTreeGuided tool;
  return tool.main(argc, argv);
}

/// @endcond

