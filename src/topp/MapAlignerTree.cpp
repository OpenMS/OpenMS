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
    /*
    // default constructor
    PeptideIdentificationsPearsonDistance(){};

    // copy constructor
    //PeptideIdentificationsPearsonDistance(const PeptideIdentificationsPearsonDistance & source);

    // destructor
    virtual ~PeptideIdentificationsPearsonDistance(){};

    // assignment operator
    PeptideIdentificationsPearsonDistance & operator=(const PeptideIdentificationsPearsonDistance & source)
    {
      if (this != &source)
      {
        DefaultParamHandler::operator=(source);
      }
      return *this;
    }
     */

    // no const SeqAndRTList because want to iterate through
    float operator()(SeqAndRTList& map_first, SeqAndRTList& map_second) const
    {
      // create vectors for both maps containing RTs of identical peptides and
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

    /*
    // calculates self similarity
    float operator()(const PeptideIdentificationsPearsonDistance & a) const
    {
      return 0;
    }

    // registers all derived products
    static void registerChildren();
    */

    static String getProductName()
    {
      return "PeptideIdentificationsPearsonDistance";
    };

  }; // end of PeptideIdentificationsPearsonDifference

  // function declarations
  template <typename MapType>
  static void loadInputMaps_(vector<MapType>& maps, StringList& ins, vector<StringList>& ms_run_paths, ConsensusMap& out_map);
  static void setUniqueIds_(vector<FeatureMap>& feature_maps);
  static void getPeptideSequences_(const vector<PeptideIdentification>& peptides, SeqAndRTList& peptide_rts, vector<double>& rts_tmp);
  static void extract_seq_and_rt_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt, vector<double>& maps_ranges);
  static void buildTree_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt, std::vector<BinaryTreeNode>& tree, vector<double>& maps_ranges);
  void treeGuidedAlignment_(const std::vector<BinaryTreeNode>& tree, vector<FeatureMap>& feature_maps,
          vector<TransformationDescription>& transformations, vector<double>& maps_ranges, ConsensusMap& out_map, String transformation_type, vector<SeqAndRTList>& maps_seqAndRT);
  static void computeTransformationsByID_(const String& transformation_type, vector<FeatureMap>& feature_maps, FeatureMap& last_map, vector<TransformationDescription>& transformations,
          const vector<Size>& trafo_order, const Param& model_params, const String& model_type);
  static void computeTransformationsByTrafo_(vector<SeqAndRTList>& maps_seqAndRT,
                                      TransformationDescription& last_trafo,
                                      vector<TransformationDescription>& transformations,
                                      const Param& model_params, const String& model_type);
  void computeConsensus_(vector<FeatureMap>& feature_maps, const vector<TransformationDescription>& transformations, ConsensusMap& out_map);
  static void storeConsensusFile_(ConsensusMap& out_map, String& out_file);
  static void storeTransformationDescriptions_(const vector<TransformationDescription>& transformations, StringList& trafos);


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
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  // getSubsectionDefaults of TOPPMapAlignerBase geht nicht, da von TOPPFeatureLinkerBase geerbt:
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
    String transformation_type = getStringOption_("transformation_type");
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    Size in_files_size = in_files.size();
    ConsensusMap out_map(in_files_size);
    vector<StringList> ms_run_paths(in_files_size);
    vector<FeatureMap> feature_maps(in_files_size);
    loadInputMaps_(feature_maps, in_files, ms_run_paths, out_map);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    setUniqueIds_(feature_maps);

    // get Peptide/ RT tuple for all features, seperated by input file
    vector<SeqAndRTList> maps_seqAndRt(in_files_size);
    // save ranges for alignment (larger rt_range -> reference)
    vector<double> maps_ranges(in_files_size);
    //extract_seq_and_rt_(feature_maps, maps_seqAndRt, maps_ranges);

    //  construct tree with pearson coefficient
    std::vector<BinaryTreeNode> tree;
    buildTree_(feature_maps, maps_seqAndRt, tree, maps_ranges);

    // print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "alignment follows tree: " << ca.newickTree(tree) << endl;

    // to store transformations
    vector<TransformationDescription> transformations(in_files_size);

    // TODO : refacture: compute transformations and consensus within treeGuidedAlignment
    treeGuidedAlignment_(tree, feature_maps, transformations, maps_ranges, out_map, transformation_type, maps_seqAndRt);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    // store transformed map
    storeConsensusFile_(out_map, out_file);

    // store transformations
    storeTransformationDescriptions_(transformations, out_trafos);

    return EXECUTION_OK;
  }
};

template<typename MapType>
void TOPPMapAlignerTree::loadInputMaps_(vector<MapType>& maps, StringList& ins, vector<StringList>& ms_run_paths,
        ConsensusMap& out_map)
        {
  FeatureXMLFile fxml_file;
  FeatureFileOptions param = fxml_file.getOptions();
  StringList ms_run_locations;

  // to save memory don't load convex hulls and subordinates
  param.setLoadSubordinates(false);
  param.setLoadConvexHull(false);
  fxml_file.setOptions(param);

  ProgressLogger progresslogger;
  progresslogger.setLogType(TOPPFeatureLinkerBase::CMD);
  progresslogger.startProgress(0, ins.size(), "loading input files");
  for (Size i = 0; i < ins.size(); ++i)
  {
    progresslogger.setProgress(i);
    fxml_file.load(ins[i], maps[i]);
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
    //OPENMS_LOG_INFO << "set " << hoho << "unique ids." << endl;
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
        TOPPMapAlignerTree::SeqAndRTList& peptide_rts, vector<double>& rts_tmp)
        {
  //vector<PeptideIdentification>::const_iterator pep_it;
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

void TOPPMapAlignerTree::extract_seq_and_rt_(const vector<FeatureMap>& feature_maps,
        vector<SeqAndRTList>& maps_seqAndRt, vector<double>& maps_ranges)
        {
  for (auto maps_it = feature_maps.begin(); maps_it != feature_maps.end(); ++maps_it)
  {
    const Size position = static_cast<Size>(distance(feature_maps.begin(), maps_it));
    double percentile10;
    double percentile90;
    vector<double> rts_tmp(maps_it->size());
    //vector<Feature>::const_iterator feature_it;
    for (auto feature_it = maps_it->begin(); maps_it->end() != feature_it; ++feature_it)
    {
      if (!feature_it->getPeptideIdentifications().empty()) {
        getPeptideSequences_(feature_it->getPeptideIdentifications(), maps_seqAndRt[position], rts_tmp);
      }
    }
    sort(rts_tmp.begin(), rts_tmp.end());

    percentile10 = rts_tmp[rts_tmp.size()*0.1];
    percentile90 = rts_tmp[rts_tmp.size()*0.9];

    maps_ranges[position] = percentile90 - percentile10;
    rts_tmp.clear();
  }
}

void TOPPMapAlignerTree::buildTree_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seqAndRt,
        std::vector<BinaryTreeNode>& tree, vector<double>& maps_ranges)
        {
  extract_seq_and_rt_(feature_maps, maps_seqAndRt, maps_ranges);
  PeptideIdentificationsPearsonDistance pepDist;
  SingleLinkage sl;
  DistanceMatrix<float> dist_matrix;
  ClusterHierarchical ch;
  ch.cluster<SeqAndRTList, PeptideIdentificationsPearsonDistance>(maps_seqAndRt, pepDist, sl, tree, dist_matrix);
}

void TOPPMapAlignerTree::treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree, vector<FeatureMap> &feature_maps,
                                              vector<TransformationDescription> &transformations,
                                              vector<double> &maps_ranges, ConsensusMap &out_map, String transformation_type,
                                              vector<SeqAndRTList>& maps_seqAndRT)
        {
  vector<TransformationDescription> trafo_tmp; // use to align
  TransformationDescription trafo_for_output;
  vector<FeatureMap> maps_transformed;
  maps_transformed = feature_maps;  // copy needed for iterations without loosing original data
  Size last_trafo = 0;  // look up transformation order in map_sets
  vector<vector<Size>> map_sets(feature_maps.size());

  Param model_params = getParam_().copy("model:", true);
  String model_type = "b_spline";
  model_params = model_params.copy(model_type+":", true);

 // vector<TransformationDescription> trafo_out(feature_maps.size());
  vector<TransformationDescription> transformations_align;
  vector<ConsensusMap> consensus_tmp(feature_maps.size() - 1);
  Size ref;
  Size to_transform;
  // helper to reorganize rt transformations
  for (Size i = 0; i < feature_maps.size(); ++i)
  {
    vector<Size> tmp;
    tmp.push_back(i);
    map_sets[i] = tmp;
  }

  MapAlignmentAlgorithmIdentification algorithm;
  Param algo_params = getParam_().copy("algorithm:", true);
  algorithm.setParameters(algo_params);
  algorithm.setLogType(log_type_);

  // perform Alignment
  for (const auto& node : tree)
  {
    vector<FeatureMap> to_align;
    //  determine the map with larger RT range for 10/90 percentile (->reference)
    if (maps_ranges[node.left_child] > maps_ranges[node.right_child]) {
        ref = node.left_child;
        to_transform = node.right_child;
        maps_ranges[node.right_child] = maps_ranges[node.left_child]; // after transformation same range for both maps
    } else {
        ref = node.right_child;
        to_transform = node.left_child;
        maps_ranges[node.left_child] = maps_ranges[node.right_child]; // after transformation same range for both maps
    }
    last_trafo = to_transform;
    // performAlignment with map as reference that has larger RT range
    to_align.push_back(maps_transformed[to_transform]);
    to_align.push_back(maps_transformed[ref]);
    // with set reference
    //algorithm.align(to_align, transformations_align, 1);

    // without set reference
    algorithm.align(to_align, transformations_align);
    /*
    if (transformations_align[0].getDataPoints()[0].first == transformations_align[0].getDataPoints()[0].second)
    {
      swap(to_transform, ref);
      swap(transformations_align[0], transformations_align[1]);
      std::cout << "nach swap: " << transformations_align[0].getDataPoints()[0].first << " : " << transformations_align[0].getDataPoints()[0].second << std::endl;
    }
    */

    // transform retention times of non-identity for next iteration
    transformations_align[0].fitModel(model_type, model_params);
    transformations_align[1].fitModel(model_type, model_params);

    // needed for following iteration steps
    MapAlignmentTransformer::transformRetentionTimes(maps_transformed[to_transform],
                                                     transformations_align[0], false);
    MapAlignmentTransformer::transformRetentionTimes(maps_transformed[ref],
                                                     transformations_align[1], false);

    // combine aligned maps, store in both, because tree always calls smaller number
    // also possible: feature_maps_transformed[smallerNumber] = ..[ref]+..[to_transform]
    // or use pointer
    maps_transformed[ref] += maps_transformed[to_transform];
    maps_transformed[ref].updateRanges();
    maps_transformed[to_transform] = maps_transformed[ref];
    trafo_for_output = transformations_align[0];

    // update transformation order for each map
    vector<Size> tmp;
    tmp = map_sets[ref];
    tmp.insert(tmp.end(), map_sets[to_transform].begin(), map_sets[to_transform].end());
    map_sets[ref] = tmp;
    map_sets[to_transform] = tmp;
    for (Size i = 0; i < feature_maps.size(); ++i)
    {
      if (i == to_transform || i == ref)
      {
        map_sets[i] = tmp;
      }
    }
    transformations_align.clear();
    to_align.clear();
  }

  // compute transformations
  if (transformation_type.empty() || transformation_type == "trafo")
  {
    computeTransformationsByTrafo_(maps_seqAndRT, trafo_for_output, transformations, model_params, model_type);
  }
  else {
    computeTransformationsByID_(transformation_type, feature_maps, maps_transformed[last_trafo], transformations, map_sets[last_trafo], model_params, model_type);
  }


  computeConsensus_(feature_maps, transformations, out_map);
}

void TOPPMapAlignerTree::computeTransformationsByTrafo_(vector<SeqAndRTList> &maps_seqAndRt,
                                                 TransformationDescription &last_trafo,
                                                 vector<TransformationDescription> &transformations,
                                                 const Param &model_params, const String &model_type) {
  ProgressLogger progresslogger;
  progresslogger.setLogType(CMD);
  progresslogger.startProgress(0, maps_seqAndRt.size(), "computing trafoXML files from trafo");
  // need to know which map was reference and which was transformed in last iteration
  Size map_id = 0;

  for (auto & map : maps_seqAndRt) {
    TransformationDescription::DataPoints trafo_data_tmp;
    auto trafoit = last_trafo.getDataPoints().begin();
    auto mapit = map.begin();
    while (trafoit != last_trafo.getDataPoints().end() && mapit != map.end()) {
      if (trafoit->note < mapit->first) {
        ++trafoit;
      } else if (trafoit->note > mapit->first) {
        ++mapit;
      } else {
        // TODO: check problems with outliers
        double rt_dist_min(std::numeric_limits<double>::infinity());
        double rt_best;
        for (auto &rt : mapit->second) {
          if (fabs(trafoit->second - rt) < rt_dist_min) {
            rt_best = rt;
          }
        }
        TransformationDescription::DataPoint point(rt_best, trafoit->second, trafoit->note);
        trafo_data_tmp.push_back(point);
        ++trafoit;
        ++mapit;
      }
    }
    transformations[map_id] = TransformationDescription(trafo_data_tmp);
    transformations[map_id].fitModel(model_type, model_params);
    ++map_id;
  }
/*
  for (Size i = 0; i < feature_maps.size(); ++i)
  {
    std::cout << "in for alles gut\n";
    auto trafo_it = last_trafos[pos[i]].getDataPoints().begin();

    for (auto feature_it = feature_maps[i].begin(); feature_it != feature_maps[i].end(); ++feature_it)
    {
      if (trafo_it != last_trafos[pos[i]].getDataPoints().end())
      {
        std::cout << "in while alles gut\n";
        TransformationDescription::DataPoints trafo_data_tmp;
        if (!feature_it->getPeptideIdentifications().empty()) {
          for (auto pep_it = feature_it->getPeptideIdentifications().begin();
               pep_it != feature_it->getPeptideIdentifications().end(); ++pep_it) {
            if (!pep_it->getHits().empty()) {
              std::cout << "hit size: " << pep_it->getHits().size() << " : " <<pep_it->getHits()[0].getScore() << std::endl;
              try {
                const string & seq = pep_it->getHits()[0].getSequence().toString();
                if (trafo_it->note < seq) {
                  ++trafo_it;
                } else if (trafo_it->note > pep_it->getHits()[0].getSequence().toString()) {
                  ++pep_it;
                } else {
                  TransformationDescription::DataPoint point(pep_it->getRT(), trafo_it->second, trafo_it->note);
                  trafo_data_tmp.push_back(point);
                  ++trafo_it;
                  ++pep_it;
                }
              }
              catch (Exception::BaseException &e)
              {
                ++pep_it;
              }
            }
          }
        }
        std::cout << "size transformations: " << transformations.size() << std::endl;
        transformations[i] = TransformationDescription(trafo_data_tmp);
        transformations[i].fitModel(model_type, model_params);
      }

    }
  }
  */
  progresslogger.endProgress();
}

void TOPPMapAlignerTree::computeTransformationsByID_(const String& transformation_type, vector<FeatureMap>& feature_maps, FeatureMap& last_map,
        vector<TransformationDescription>& transformations, const vector<Size>& trafo_order, const Param& model_params,
        const String& model_type) {
  ProgressLogger progresslogger;
  progresslogger.setLogType(CMD);
  progresslogger.startProgress(0, 1, "computing trafoXML files by id");
  auto last_map_it = last_map.begin();
  for (auto & map_idx : trafo_order)
  {
    TransformationDescription::DataPoints trafo_data_tmp;
    for (auto & features : feature_maps[map_idx])
    {
      if (features.getUniqueId() == last_map_it->getUniqueId())
      {
        if (transformation_type == "peptides")
        {
          auto last_map_pep_it = last_map_it->getPeptideIdentifications().begin();
          auto pep_it = features.getPeptideIdentifications().begin();
          while (last_map_pep_it != last_map_it->getPeptideIdentifications().end() &&
                 pep_it != features.getPeptideIdentifications().end())
          {
            if (last_map_pep_it->getHits()[0].getSequence() == pep_it->getHits()[0].getSequence())
            {
              TransformationDescription::DataPoint point(pep_it->getRT(), last_map_it->getRT(), pep_it->getHits()[0].getSequence().toString());
              trafo_data_tmp.push_back(point);
            }
            else{
              OPENMS_LOG_INFO << "peptide identification hits don't have the same sequence" << std::endl;
            }
            ++last_map_pep_it;
            ++pep_it;
          }
        }
        else if (transformation_type == "features")
        {
          TransformationDescription::DataPoint point(features.getRT(), last_map_it->getRT(), features.getUniqueId());
          trafo_data_tmp.push_back(point);
        }


      } else{
        OPENMS_LOG_INFO << "features to compare don't have the same unique id" << endl;
      }
      ++last_map_it;
    }
    transformations[map_idx] = TransformationDescription(trafo_data_tmp);
    transformations[map_idx].fitModel(model_type, model_params);
    trafo_data_tmp.clear();
  }
  progresslogger.endProgress();
}

void TOPPMapAlignerTree::computeConsensus_(vector<FeatureMap>& feature_maps, const vector<TransformationDescription>& transformations, ConsensusMap& out_map) {

  ProgressLogger progresslogger;
  progresslogger.setLogType(CMD);
  progresslogger.startProgress(0, 1, "computing consensus map");
  auto trafo = transformations.begin();
  for (auto& map : feature_maps)
  {
    MapAlignmentTransformer::transformRetentionTimes(map,
                                                     *trafo, false);
    map.updateRanges(); // without: LeakSanitizer detects memory leaks
    ++trafo;
  }
  FeatureGroupingAlgorithmKD link_feature_maps;
  Param p = link_feature_maps.getDefaults();
  p.setValue("warp:enabled", "true"); // no additional rt transformation by FeatureLinker
  link_feature_maps.setParameters(p);
  link_feature_maps.group(feature_maps, out_map);

  //assign unique ids
  out_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

  // annotate output with data processing info
  addDataProcessing_(out_map,
                     getProcessingInfo_(DataProcessing::FEATURE_GROUPING));


  // sort list of peptide identifications in each consensus feature by map index
  out_map.sortPeptideIdentificationsByMapIndex();
  progresslogger.endProgress();
}

void TOPPMapAlignerTree::storeConsensusFile_(ConsensusMap& out_map, String& out_file) {
  ConsensusXMLFile cxml_file;

  ProgressLogger progresslogger;
  progresslogger.setLogType(CMD);
  progresslogger.startProgress(0, 1, "writing output file");
  // annotate output with data processing info:
  //addDataProcessing_(maps[trafo_order.back()],
  //                   getProcessingInfo_(DataProcessing::ALIGNMENT));
  //maps[trafo_order.back()].setUniqueId();
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

void TOPPMapAlignerTree::storeTransformationDescriptions_(const vector<TransformationDescription>& transformations,
                                                          StringList& trafos) {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(CMD);
    progresslogger.startProgress(0, trafos.size(),
                                 "writing transformation files");
    for (Size i = 0; i < trafos.size(); ++i)
    {
        TransformationXMLFile().store(trafos[i], transformations[i]);
    }
    progresslogger.endProgress();
}


int main(int argc, const char** argv)
{
  TOPPMapAlignerTree tool;
  return tool.main(argc, argv);
}

/// @endcond
