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
  static void loadInputMaps_(vector<MapType>& maps, const StringList& in_files, vector<StringList>& ms_run_locations);
  static void setUniqueIds_(vector<FeatureMap>& feature_maps);
  static void getPeptideSequences_(const vector<PeptideIdentification>& peptides, SeqAndRTList& peptide_rts);
  static void extract_seq_and_rt_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt);
  static void buildTree_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt, std::vector<BinaryTreeNode>& tree);
  void treeGuidedAlignment_(const std::vector<BinaryTreeNode>& tree, vector<ConsensusMap>& consensus_maps,
                            TransformationDescription& transformation, ConsensusMap& out_map, vector<FeatureMap>& feature_maps, Size& last_map_idx);
  void createTransformationFiles_(const vector<FeatureMap> &feature_maps, const ConsensusMap &out_map,
                                  const TransformationDescription &last_trafo_descr,
                                  vector<TransformationDescription> &transformations);
  void storeConsensusFile_(ConsensusMap& out_map, String& out_file);


  void registerOptionsAndFlags_() override
  {
    TOPPFeatureLinkerBase::registerOptionsAndFlags_();
    //registerInputFileList_("in", "<files>", ListUtils::create<String>(""), "Input files", true);
    //setValidFormats_("in", ListUtils::create<String>("featureXML"));
    //registerOutputFile_("out", "<file>", "", "Output file", true);
    //setValidFormats_("out", ListUtils::create<String>("consensusXML"));
    registerOutputFileList_("trafo_out", "<files>", StringList(), "Transformation output files. This option or 'out' has to be provided; they can be used together.", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    registerStringOption_("transformation_type", "string", "trafo", "Option to decide transformation path during alignment.", false);
    setValidStrings_("transformation_type",ListUtils::create<String>("trafo,features,peptides"));
    registerStringOption_("fl_rt_transform", "string", "true", "With true the FeatureLinkerUnlabeldKD transforms retention times of input files.", false);
    setValidStrings_("fl_rt_transform", ListUtils::create<String>("true,false"));
    //registerStringOption_("keep_subelements", "string", "true", "kepp subelements after grouping features", false);
    //setValidStrings_("keep_subelements", ListUtils::create<String>("true,false"));
    registerSubsection_("align_algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
    registerSubsection_("linker_algorithm", "FeatureGroupingAlgorithm");
  }

  Param getSubsectionDefaults_(const String&  section) const override
  {
    if (section == "align_algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return TOPPMapAlignerBase::getModelDefaults("b_spline");
    }
    if (section == "linker_algorithm")
    {
      FeatureGroupingAlgorithmKD algo;
      Param p = algo.getParameters();
      return p;
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
    vector<StringList> ms_run_paths(in_files_size);
    vector<FeatureMap> feature_maps(in_files_size);
    loadInputMaps_(feature_maps, in_files, ms_run_paths);
    std::cout << "in_files size: " << in_files.size() << " maps size: " << feature_maps.size() << std::endl;

    setUniqueIds_(feature_maps);

    // ------------- convert to ConsensusMap----------------------
    vector<ConsensusMap> consensus_maps(feature_maps.size());
    Int max_num_peaks_considered_ = 500; // convert uses size of map, if this value is higher
    for (Size i = 0; i < feature_maps.size(); ++i)
    {
      // unique ids and ranges updated by convert
      MapConversion::convert(0, feature_maps[i], consensus_maps[i], max_num_peaks_considered_);
      consensus_maps[i].getColumnHeaders()[0].filename = ms_run_paths[i].front();
      consensus_maps[i].getColumnHeaders()[0].unique_id = feature_maps[i].getUniqueId();
      //consensus_maps[i].updateUniqueIdToIndex();

      std::cout << "map headers: " << consensus_maps[i].getColumnHeaders()[0].filename << std::endl;
      std::cout << "map headers: " << consensus_maps[i].getColumnHeaders()[0].size << std::endl;
      std::cout << "map headers: " << consensus_maps[i].getColumnHeaders()[0].unique_id << std::endl;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // get Peptide/ RT tuple for all features, seperated by input file
    vector<SeqAndRTList> maps_seqAndRt(in_files_size);

    //  construct tree with pearson coefficient
    std::vector<BinaryTreeNode> tree;
    buildTree_(feature_maps, maps_seqAndRt, tree);

    // print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "alignment follows tree: " << ca.newickTree(tree) << endl;

    // to store last transformation
    TransformationDescription transformation;


    // TODO : refacture: compute transformations and consensus within treeGuidedAlignment
    Size last_map_idx;
    treeGuidedAlignment_(tree, consensus_maps, transformation, out_map, feature_maps, last_map_idx);

    vector<TransformationDescription> transformations;
    createTransformationFiles_(feature_maps, out_map, transformation, transformations);
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    // store transformed map
    /*
    std::cout << "UID header 0: " << out_map.getColumnHeaders()[0].unique_id << std::endl;
    out_map.getColumnHeaders()[0].unique_id = feature_maps[0].getUniqueId();
    out_map.getColumnHeaders()[0].size = feature_maps[0].size();
    out_map.getColumnHeaders()[0].filename = ms_run_paths[0][0];

    out_map.getColumnHeaders()[1].unique_id = feature_maps[1].getUniqueId();
    out_map.getColumnHeaders()[1].size = feature_maps[1].size();
    out_map.getColumnHeaders()[1].filename = ms_run_paths[1][0];
     */

    storeConsensusFile_(consensus_maps[last_map_idx], out_file);

    // store transformations
    //storeTransformationDescriptions_(transformations, out_trafos);

    return EXECUTION_OK;
  }
};

template<typename MapType>
void TOPPMapAlignerTree::loadInputMaps_(vector<MapType>& maps, const StringList& in_files, vector<StringList>& ms_run_locations)
{
  FeatureXMLFile fxml_file;
  FeatureFileOptions param = fxml_file.getOptions();

  // to save memory don't load convex hulls and subordinates
  param.setLoadSubordinates(false);
  param.setLoadConvexHull(false);
  fxml_file.setOptions(param);

  ProgressLogger progresslogger;
  Size progress = 0;
  progresslogger.setLogType(TOPPFeatureLinkerBase::CMD);
  progresslogger.startProgress(0, in_files.size(), "loading input files");
  for (Size i = 0; i < in_files.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.load(in_files[i], maps[i]);
    maps[i].getPrimaryMSRunPath(ms_run_locations[i]);

    // associate mzML file with map i in consensusXML
    if (ms_run_locations[i].size() > 1 || ms_run_locations[i].empty())
    {
      OPENMS_LOG_WARN << "Exactly one MS runs should be associated with a FeatureMap. "
                      << ms_run_locations[i].size()
                      << " provided." << endl;
    }

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
    progresslogger.setProgress(progress++);
  }
  progresslogger.endProgress();
}

void TOPPMapAlignerTree::setUniqueIds_(vector<FeatureMap> &feature_maps) {
  FeatureMap maps;
  for (auto& map : feature_maps)
  {
    maps += map;
  }
  Size setUID = maps.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  Size resolveUID = maps.resolveUniqueIdConflicts();
  std::cout << "setUID: " << setUID << " resolve: " << resolveUID << std::endl;
  //if (maps.applyMemberFunction(&UniqueIdInterface::setUniqueId))
  //{
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
  //}
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

void TOPPMapAlignerTree::treeGuidedAlignment_(const std::vector<BinaryTreeNode>& tree, vector<ConsensusMap>& consensus_maps,
        TransformationDescription& transformation, ConsensusMap& out_map, vector<FeatureMap>& feature_maps, Size& last_map_idx)
        {
  MapAlignmentAlgorithmIdentification algorithm;
  algorithm.setLogType(CMD);

  Param model_params = getParam_().copy("model:", true);
  String model_type = "b_spline";
  model_params = model_params.copy(model_type+":", true);

  FeatureGroupingAlgorithmKD link_feature_maps;
  Param p = getParam_().copy("linker_algorithm:", true);
  //hard coded, weil ich es anders nicht hin bekomme :(
  p.setValue("keep_subelements", true);
  p.setValue("nr_partitions", 1);
  p.setValue("warp:enabled", "false");
  p.setValue("mz_unit", "Da");
  p.setValue("warp:mz_tol", 0.3);
  p.setValue("link:rt_tol", 100.0);
  p.setValue("link:mz_tol", 0.3);

  link_feature_maps.setParameters(p);

  //last_map_idx = 0; // use to know where last consensus is stored
  // align maps tree guided
  vector<vector<Size>> order(feature_maps.size());
  for (Size i = 0; i < feature_maps.size(); ++i)
  {
    order[i].push_back(i);
  }
  for (const auto & node : tree)
  {
    vector<ConsensusMap> to_align;
    vector<TransformationDescription> tmp_trafoDesc;
    order[node.left_child].insert(order[node.left_child].end(), order[node.right_child].begin(), order[node.right_child].end());
    to_align.push_back(consensus_maps[node.left_child]);
    to_align.push_back(consensus_maps[node.right_child]);
    algorithm.align(to_align, tmp_trafoDesc);
    tmp_trafoDesc[0].fitModel(model_type, model_params);
    tmp_trafoDesc[1].fitModel(model_type, model_params);
    for (Size i = 0; i < to_align.size(); ++i)
    {

      MapAlignmentTransformer::transformRetentionTimes(to_align[i], tmp_trafoDesc[i]);
      to_align[i].updateRanges();

    }
    // use featureGrouping to get consensus of alignment
    ConsensusMap cons_tmp;
    link_feature_maps.group(to_align, cons_tmp);
    consensus_maps[node.left_child].clear();
    consensus_maps[node.left_child] = cons_tmp;
    /*
    for (auto & map_idx : order[node.left_child])
    {
      StringList ms_runs;
      feature_maps[map_idx].getPrimaryMSRunPath(ms_runs);
      consensus_maps[node.left_child].getColumnHeaders()[map_idx].filename =ms_runs[0];
      consensus_maps[node.left_child].getColumnHeaders()[map_idx].size = feature_maps[map_idx].size();
      consensus_maps[node.left_child].getColumnHeaders()[map_idx].unique_id = feature_maps[map_idx].getUniqueId();
    }
    Size ensuredIDs = consensus_maps[node.left_child].applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
    std::cout << "ensured nach grouping: " << ensuredIDs << std::endl;
     */
    last_map_idx = node.left_child;  // need to know position of last consensus
    //consensus_maps[node.right_child] = consensus_maps[node.left_child]; // not needed because tree always links to left child?
    link_feature_maps.transferSubelements(to_align, consensus_maps[node.left_child]);
    consensus_maps[node.left_child].sortPeptideIdentificationsByMapIndex();
  }

  out_map = consensus_maps[last_map_idx];
  std::cout << "resolveUniquIDs in out map: " << out_map.resolveUniqueIdConflicts() << std::endl;
}

void TOPPMapAlignerTree::createTransformationFiles_(const vector<FeatureMap> &feature_maps, const ConsensusMap &out_map,
                                                    const TransformationDescription &last_trafo_descr,
                                                    vector<TransformationDescription> &transformations) {
  for (const auto & map : feature_maps)
  {
    UInt64 unique_id_map = map.getUniqueId();
    for (auto & header : out_map.getColumnHeaders())
    {
      //if (header.unique_id)
    }
  }

}

void TOPPMapAlignerTree::storeConsensusFile_(ConsensusMap& out_map, String& out_file) {
  ConsensusXMLFile cxml_file;

  ProgressLogger progresslogger;
  progresslogger.setLogType(CMD);
  progresslogger.startProgress(0, 1, "writing output file");
  // annotate output with data processing info:
  //addDataProcessing_(maps[trafo_order.back()],
  //                   getProcessingInfo_(DataProcessing::ALIGNMENT));
  //assign unique ids
  //Size ensuredIDs = out_map.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  // annotate output with data processing info

  //out_map.sortPeptideIdentificationsByMapIndex();
  addDataProcessing_(out_map,
                     getProcessingInfo_(DataProcessing::FEATURE_GROUPING));
  cxml_file.store(out_file, out_map);
  progresslogger.endProgress();

  // some statistics
  map<Size, UInt> num_consfeat_of_size;
  for (ConsensusMap::const_iterator cmit = out_map.begin();
       cmit != out_map.end(); ++cmit)
  {
    ++num_consfeat_of_size[cmit->size()];
  }

  OPENMS_LOG_INFO << "Number of consensus features:" << endl;
  for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin();
       i != num_consfeat_of_size.rend(); ++i)
  {
    OPENMS_LOG_INFO << "  of size " << setw(2) << i->first << ": " << setw(6)
                    << i->second << endl;
  }
  OPENMS_LOG_INFO << "  total:      " << setw(6) << out_map.size() << endl;
}


int main(int argc, const char** argv)
{
  TOPPMapAlignerTree tool;
  return tool.main(argc, argv);
}

/// @endcond