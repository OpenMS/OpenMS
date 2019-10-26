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

#include "FeatureLinkerBase.cpp"
// calculate pearson distance
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
// create binary tree
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
// align maps and genereate output
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
// store transformations
#include <OpenMS/FORMAT/TransformationXMLFile.h>

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

class TOPPMapAlignerTree :
        public TOPPFeatureLinkerBase
{

public:
  TOPPMapAlignerTree() :
          TOPPFeatureLinkerBase("MapAlignerTree", "Tree guided correction of retention time distortions between maps.")
  {
  }

private:
  /// Type to store retention times given for individual peptide sequence
  typedef std::map<String, DoubleList> SeqAndRTList;

  class PeptideIdentificationsPearsonDistance
  {
  public:
    // no const SeqAndRTList because want to iterate through
    float operator()(SeqAndRTList& map_first, SeqAndRTList& map_second) const
    {
      // create vectors for both maps containing RTs of peptides with same sequence and
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
          //if (pep1_it->second.size()+pep2_it->second.size() >2) std::cout << pep1_it->second.front() << " " << pep1_it->second.back() << " " << pep2_it->second.front()<< " " <<pep2_it->second.back() << " " << std:: endl;
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

    static String getProductName()
    {
      return "PeptideIdentificationsPearsonDistance";
    };

  }; // end of PeptideIdentificationsPearsonDifference

  // function declarations
  template <typename MapType>
  static void loadInputMaps_(vector<MapType>& maps, const StringList& in_files, ConsensusMap& out_map);
  static void setUniqueIds_(vector<FeatureMap>& feature_maps);
  static void getPeptideSequences_(const vector<PeptideIdentification>& peptides, SeqAndRTList& peptide_rts);
  static void extract_seq_and_rt_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt);
  static void buildTree_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt, std::vector<BinaryTreeNode>& tree);
  void treeGuidedAlignment_(const std::vector<BinaryTreeNode>& tree, const vector<ConsensusMap>& consensus_maps,
                            vector<TransformationDescription>& transformations, ConsensusMap& out_map);


  void registerOptionsAndFlags_() override
  {
    //TOPPFeatureLinkerBase::registerOptionsAndFlags_();
    registerInputFileList_("in", "<files>", ListUtils::create<String>(""), "Input files", true);
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out", "<file>", "", "Output file", true);
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFileList_("trafo_out", "<files>", StringList(), "Transformation output files. This option or 'out' has to be provided; they can be used together.", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    registerStringOption_("transformation_type", "string", "trafo", "Option to decide transformation path during alignment.", false);
    setValidStrings_("transformation_type",ListUtils::create<String>("trafo,features,peptides"));
    registerStringOption_("fl_rt_transform", "string", "true", "With true the FeatureLinkerUnlabeldKD transforms retention times of input files.", false);
    setValidStrings_("fl_rt_transform", ListUtils::create<String>("true,false"));
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  Param getSubsectionDefaults_(const String&  section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return TOPPMapAlignerBase::getModelDefaults("b_spline");
    }

    return Param(); // this shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    //ExitCodes ret = checkParameters_();
    //if (ret != EXECUTION_OK) return ret;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in_files = getStringList_("in");
    String out_file = getStringOption_("out");
    StringList out_trafos = getStringList_("trafo_out");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    Size in_files_size = in_files.size();
    ConsensusMap out_map(in_files_size);
    vector<FeatureMap> feature_maps(in_files_size);
    loadInputMaps_(feature_maps, in_files, out_map);

    // ------------- convert to ConsensusMap----------------------
    vector<ConsensusMap> consensus_maps(feature_maps.size());
    Int max_num_peaks_considered_ = 10000; // convert uses size of map, if this value is higher
    for (Size i = 0; i < feature_maps.size(); ++i)
    {
      // unique ids kept and ranges updated by convert
      MapConversion::convert(1, feature_maps[i], consensus_maps[i], max_num_peaks_considered_);
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    setUniqueIds_(feature_maps);

    // get Peptide/ RT tuple for all features, seperated by input file
    vector<SeqAndRTList> maps_seqAndRt(in_files_size);

    //  construct tree with pearson coefficient
    std::vector<BinaryTreeNode> tree;
    buildTree_(feature_maps, maps_seqAndRt, tree);

    // print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "alignment follows tree: " << ca.newickTree(tree) << endl;

    // to store transformations
    vector<TransformationDescription> transformations(in_files_size);


    // TODO : refacture: compute transformations and consensus within treeGuidedAlignment
    treeGuidedAlignment_(tree, consensus_maps, transformations, out_map);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    // store transformed map
    //storeConsensusFile_(out_map, out_file);

    // store transformations
    //storeTransformationDescriptions_(transformations, out_trafos);

    return EXECUTION_OK;
  }
};

template<typename MapType>
void TOPPMapAlignerTree::loadInputMaps_(vector<MapType>& maps, const StringList& in_files,
                                        ConsensusMap& out_map)
{
  vector<StringList> ms_run_paths(in_files.size());
  FeatureXMLFile fxml_file;
  FeatureFileOptions param = fxml_file.getOptions();
  StringList ms_run_locations;

  // to save memory don't load convex hulls and subordinates
  param.setLoadSubordinates(false);
  param.setLoadConvexHull(false);
  fxml_file.setOptions(param);

  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPFeatureLinkerBase::CMD);
  progresslogger.startProgress(0, in_files.size(), "loading input files");
  for (Size i = 0; i < in_files.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.load(in_files[i], maps[i]);
    maps[i].getPrimaryMSRunPath(ms_run_paths[i]);
    if (ms_run_paths[i].size() > 1 || ms_run_paths[i].empty())
    {
      OPENMS_LOG_WARN << "Exactly one MS runs should be associated with a FeatureMap. "
                      << ms_run_paths[i].size()
                      << " provided." << endl;
    }
    else
    {
      out_map.getColumnHeaders()[i].filename = ms_run_paths[i].front();
    }
    out_map.getColumnHeaders()[i].size = maps[i].size();
    out_map.getColumnHeaders()[i].unique_id = maps[i].getUniqueId();

    // copy over information on the primary MS run
    //ms_run_locations.insert(ms_run_locations[i].end(), ms_run_paths.begin(), ms_run_paths.end());

    // to save memory, remove convex hulls, subordinates:
    for (auto it = maps[i].begin(); it != maps[i].end();
         ++it)
    {
      String adduct;
      //exception: addduct information
      if (it->metaValueExists("dc_charge_adducts"))
      {
        adduct = it->getMetaValue("dc_charge_adducts");
      }
      it->getSubordinates().clear();
      it->getConvexHulls().clear();
      it->clearMetaInfo();
      if (!adduct.empty())
      {
        it->setMetaValue("dc_charge_adducts", adduct);
      }
    }
    maps[i].updateRanges();
  }
  progresslogger.endProgress();
}

void TOPPMapAlignerTree::setUniqueIds_(vector<FeatureMap> &feature_maps) {
  FeatureMap maps;
  for (auto& map : feature_maps)
  {
    maps += map;
  }
  if (maps.applyMemberFunction(&UniqueIdInterface::setUniqueId))
  {
    auto maps_it = maps.begin();
    for (auto& map : feature_maps)
    {
      for (auto& feature : map)
      {
        if (feature.getUniqueId() != maps_it->getUniqueId())
        {
          feature.setUniqueId(maps_it->getUniqueId());
        }
        ++maps_it;
      }
    }
  }
}

void TOPPMapAlignerTree::getPeptideSequences_(const vector<PeptideIdentification>& peptides,
                                              TOPPMapAlignerTree::SeqAndRTList& peptide_rts)
{
  for (const auto & peptide : peptides)
  {
    if (!peptide.getHits().empty())
    {
      const String& sequence = peptide.getHits()[0].getSequence().toString();
      double rt = peptide.getRT();
      peptide_rts[sequence].push_back(rt);
    }
  }
}

void TOPPMapAlignerTree::extract_seq_and_rt_(const vector<FeatureMap>& feature_maps,
                                             vector<SeqAndRTList>& maps_seqAndRt)
{
  for (auto maps_it = feature_maps.begin(); maps_it != feature_maps.end(); ++maps_it)
  {
    const Size position = static_cast<Size>(distance(feature_maps.begin(), maps_it));
    for (auto feature_it = maps_it->begin(); maps_it->end() != feature_it; ++feature_it)
    {
      if (!feature_it->getPeptideIdentifications().empty()) {
        getPeptideSequences_(feature_it->getPeptideIdentifications(), maps_seqAndRt[position]);
      }
    }
  }
}

void TOPPMapAlignerTree::buildTree_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt,
                                    std::vector<BinaryTreeNode>& tree)
{
  extract_seq_and_rt_(feature_maps, maps_seqAndRt);
  PeptideIdentificationsPearsonDistance pepDist;
  SingleLinkage sl;
  DistanceMatrix<float> dist_matrix;
  ClusterHierarchical ch;
  ch.cluster<SeqAndRTList, PeptideIdentificationsPearsonDistance>(maps_seqAndRt, pepDist, sl, tree, dist_matrix);
}

void TOPPMapAlignerTree::treeGuidedAlignment_(const std::vector<BinaryTreeNode>& tree, const vector<ConsensusMap>& consensus_maps,
        vector<TransformationDescription>& transformations, ConsensusMap& out_map)
        {
  MapAlignmentAlgorithmIdentification algorithm;
  Param algo_params = getParam_().copy("algorithm:", true);
  algorithm.setParameters(algo_params);
  algorithm.setLogType(CMD);

  Param model_params = getParam_().copy("model:", true);
  String model_type = "b_spline";
  model_params = model_params.copy(model_type+":", true);

  FeatureGroupingAlgorithmKD link_feature_maps;
  Param p = link_feature_maps.getDefaults();
  p.setValue("warp:enabled", true); // no additional rt transformation by FeatureLinker
  link_feature_maps.setParameters(p);

  Size last_node; // use to know where last consensus is stored
  // align maps tree guided
  for (const auto & node : tree)
  {
    vector<ConsensusMap> to_align;
    vector<TransformationDescription> tmp_trafoDesc;
    to_align.push_back(consensus_maps[node.left_child]);
    to_align.push_back(consensus_maps[node.right_child]);
    algorithm.align(to_align, tmp_trafoDesc);
    tmp_trafoDesc[0].fitModel(model_type, model_params);
    tmp_trafoDesc[1].fitModel(model_type, model_params);
    for (Size i = 0; i < to_align.size(), ++i)
    {
      MapAlignmentTransformer::transformRetentionTimes(to_align[i], tmp_trafoDesc[i]);
      to_align[i].updateRanges();
    }
    // use featureGrouping to get consensus od alignment
    link_feature_maps.group(to_align, consensus_maps[node.left_child]);
    last_node = node.left_child;  // need to know position of last consensus
    //consensus_maps[node.right_child] = consensus_maps[node.left_child]; // not needed because tree always links to left child?
  }

  out_map = consensus_maps[last_node];
}


int main(int argc, const char** argv)
{
  TOPPMapAlignerTree tool;
  return tool.main(argc, argv);
}

/// @endcond