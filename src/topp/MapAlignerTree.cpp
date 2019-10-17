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
  {}

private:
  /// Type to store retention times given for individual peptide sequence
  typedef std::map<String, DoubleList> SeqAndRTList;

  class PeptideIdentificationsPearsonDistance
  {
  public:
    /*
    // default constructor
    PeptideIdentificationsPearsonDistance();

    // copy constructor
    PeptideIdentificationsPearsonDistance(const PeptideIdentificationsPearsonDistance & source);

    // destructor
    ~PeptideIdentificationsPearsonDistance() = default_delete;

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

    float operator()(SeqAndRTList& map_first, SeqAndRTList& map_second) const
    {
      // create vectors for both maps containing RTs of identical peptides and
      // get union and intercept amount of peptides
      SeqAndRTList::iterator pep1_it = map_first.begin();
      SeqAndRTList::iterator pep2_it = map_second.begin();
      vector<double> intercept_rts1;
      vector<double> intercept_rts2;
      int union_size = 0;
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

    static String getProductName();

  }; // end of PeptideIdentificationsPearsonDifference

  // function declarations
  template <typename MapType, typename FileType>
  void loadInputMaps_(vector<MapType> &maps, StringList &ins, FileType &input_file);
  static void getPeptideSequences_(vector<PeptideIdentification> &peptides, SeqAndRTList &peptide_rts, vector<double> &rts_tmp);
  static void extract_seq_and_rt_(vector<FeatureMap> &feature_maps, vector<SeqAndRTList> &maps_seqAndRt, vector<double> &maps_ranges);
  static void buildTree_(vector<SeqAndRTList> &maps_seqAndRt, std::vector<BinaryTreeNode> &tree);
  void treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree, vector<FeatureMap> &feature_maps_transformed, vector<TransformationDescription> &transformations_tmp, vector<double> &maps_ranges, vector<int> &trafo_order, const Param &model_params, const String &model_type);
  static void computeTransformations_(vector<SeqAndRTList> &maps_seqAndRt, const vector<TransformationDescription> &trafo_tmp, vector<TransformationDescription> &transformations, const Param &model_params, const String &model_type);

  void storeConsensusFile_(vector<FeatureMap> &feature_maps, vector<int> &trafo_order, FeatureXMLFile &fxml_file, String out_file)
  {
      OPENMS_LOG_INFO << "writing output file" << endl;
      // annotate output with data processing info:
      addDataProcessing_(feature_maps[trafo_order.back()],
                         getProcessingInfo_(DataProcessing::ALIGNMENT));
      feature_maps[trafo_order.back()].setUniqueId();
      fxml_file.store(out_file, feature_maps[trafo_order.back()]);
  }

  void storeTransformationDescriptions_(const vector<TransformationDescription>&
                                        transformations, StringList& trafos)
  {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, trafos.size(),
                                 "writing transformation files");
    for (Size i = 0; i < transformations.size(); ++i)
    {
      TransformationXMLFile().store(trafos[i], transformations[i]);
    }
    progresslogger.endProgress();
  }

  void registerOptionsAndFlags_() override
  {
    String formats = "featureXML";
    //TOPPFeatureLinkerBase::registerOptionsAndFlags_();
    registerInputFileList_("in", "<files>", ListUtils::create<String>(""), "Input files", true);
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out", "<file>", "", "Output file", true);
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFileList_("trafo_out", "<files>", StringList(), "Transformation output files. This option or 'out' has to be provided; they can be used together.", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  // getSubsectionDefaults of TOPPMapAlignerBase geht nicht, da von TOPPFeatureLinkerBase geerbt:
  Param getSubsectionDefaults_(const String & section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      // deswegen ist keine Definition in der INI noetig.. somit aber auch nicht parametrisierbar!?
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
    vector<FeatureMap> feature_maps(in_files_size);
    FeatureXMLFile fxml_file;
    //fxml_file.getOptions();
    loadInputMaps_(feature_maps, in_files, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // get Peptide/ RT tuple for all features, seperated by input file
    vector<SeqAndRTList> maps_seqAndRt(in_files_size);
    // save ranges for alignment (larger rt_range -> reference)
    vector<double> maps_ranges(in_files_size);
    extract_seq_and_rt_(feature_maps, maps_seqAndRt, maps_ranges);

    //  construct tree with pearson coefficient
    std::vector<BinaryTreeNode> tree;
    buildTree_(maps_seqAndRt, tree);

    // to print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "alignment follows tree: " << ca.newickTree(tree) << endl;

    // to store transformations
    vector<TransformationDescription> trafo_tmp;
    vector<FeatureMap> feature_maps_transformed = feature_maps; //copy needed for computation of transformations
    vector<int> trafo_order;

    // perform Alignment
    Param model_params = getParam_().copy("model:", true);
    String model_type = "b_spline";
    model_params = model_params.copy(model_type+":", true);

    treeGuidedAlignment_(tree, feature_maps_transformed, trafo_tmp, maps_ranges, trafo_order, model_params, model_type);

    // generate transformations for each map
    vector<TransformationDescription> transformations;
    computeTransformations_(maps_seqAndRt, trafo_tmp, transformations, model_params, model_type);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    // store transformed map
    storeConsensusFile_(feature_maps_transformed, trafo_order, fxml_file, out_file);

    // store transformations
    storeTransformationDescriptions_(transformations, out_trafos);

    return EXECUTION_OK;
  }

};

String TOPPMapAlignerTree::PeptideIdentificationsPearsonDistance::getProductName() {
    return "PeptideIdentificationsPearsonDistance";
}

template<typename MapType, typename FileType>
void TOPPMapAlignerTree::loadInputMaps_(vector<MapType> &maps, StringList &ins, FileType &input_file) {
    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPFeatureLinkerBase::log_type_);
    progresslogger.startProgress(0, ins.size(), "loading input files");
    for (Size i = 0; i < ins.size(); ++i)
    {
        progresslogger.setProgress(i);
        input_file.load(ins[i], maps[i]);
    }
    progresslogger.endProgress();
}

void TOPPMapAlignerTree::getPeptideSequences_(vector<PeptideIdentification> &peptides,
                                              TOPPMapAlignerTree::SeqAndRTList &peptide_rts, vector<double> &rts_tmp) {
    vector<PeptideIdentification>::iterator pep_it;
    for (pep_it = peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
        if (!pep_it->getHits().empty())
        {
            const String& sequence = pep_it->getHits()[0].getSequence().toString();
            double rt = pep_it->getRT();
            peptide_rts[sequence].push_back(rt);
            rts_tmp.push_back(rt);
        }
    }
}

void TOPPMapAlignerTree::extract_seq_and_rt_(vector<FeatureMap> &feature_maps, vector<SeqAndRTList> &maps_seqAndRt, vector<double> &maps_ranges) {
    for (auto maps_it = feature_maps.begin(); maps_it != feature_maps.end(); ++maps_it)
    {
        Size position = static_cast<Size>(distance(feature_maps.begin(), maps_it));
        double percentile10;
        double percentile90;
        vector<double> rts_tmp(maps_it->size());
        vector<Feature>::iterator feature_it;
        for (feature_it = maps_it->begin();
             maps_it->end() != feature_it; ++feature_it)
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

void TOPPMapAlignerTree::buildTree_(vector<SeqAndRTList> &maps_seqAndRt, std::vector<BinaryTreeNode> &tree) {
    PeptideIdentificationsPearsonDistance pepDist;
    SingleLinkage sl;
    DistanceMatrix<float> dist_matrix;
    ClusterHierarchical ch;
    ch.cluster<SeqAndRTList, PeptideIdentificationsPearsonDistance>(maps_seqAndRt, pepDist, sl, tree, dist_matrix);
}

void TOPPMapAlignerTree::treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree,
                                              vector<FeatureMap> &feature_maps_transformed,
                                              vector<TransformationDescription> &transformations_tmp,
                                              vector<double> &maps_ranges, vector<int> &trafo_order,
                                              const Param &model_params, const String &model_type) {
    vector<TransformationDescription> transformations_align;
    vector<FeatureMap> to_align;
    Size ref;
    Size to_transform;

    MapAlignmentAlgorithmIdentification algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    for (BinaryTreeNode node : tree)
    {
        //  determine the map with larger RT range for 10/90 percentile (->reference)
        if (maps_ranges[node.left_child] > maps_ranges[node.right_child]) {
            ref = node.left_child;
            to_transform = node.right_child;
            trafo_order.push_back(to_transform);
            maps_ranges[node.right_child] = maps_ranges[node.left_child]; // after transformation same range for both maps
        } else {
            ref = node.right_child;
            to_transform = node.left_child;
            trafo_order.push_back(to_transform);
            maps_ranges[node.left_child] = maps_ranges[node.right_child]; // after transformation same range for both maps
        }
        // performAlignment with map as reference that has larger RT range
        algorithm.setReference(feature_maps_transformed[ref]);
        to_align.push_back(feature_maps_transformed[to_transform]);
        to_align.push_back(feature_maps_transformed[ref]);
        algorithm.align(to_align, transformations_align, 1);

        // transform retention times of non-identity for next iteration
        transformations_align[0].fitModel(model_type, model_params);
        // needed for following iteration steps
        MapAlignmentTransformer::transformRetentionTimes(feature_maps_transformed[to_transform],
                transformations_align[0], true);

        // combine aligned maps, store in both, because tree always calls smaller number
        // also possible: feature_maps_transformed[smallerNumber] = ..[ref]+..[to_transform] or use pointer
        feature_maps_transformed[ref] += feature_maps_transformed[to_transform];
        feature_maps_transformed[ref].updateRanges();
        //feature_maps_transformed[ref].updateUniqueIdToIndex();
        feature_maps_transformed[to_transform] = feature_maps_transformed[ref];
        //feature_maps_transformed[ref].updateUniqueIdToIndex();

        // alternative for insertion using +=:


        transformations_tmp.push_back(transformations_align[0]);

        transformations_align.clear();
        to_align.clear();
    }
}

void TOPPMapAlignerTree::computeTransformations_(vector<SeqAndRTList> &maps_seqAndRt,
                                                 const vector<TransformationDescription> &trafo_tmp,
                                                 vector<TransformationDescription> &transformations,
                                                 const Param &model_params, const String &model_type) {
    for (auto & mapsit : maps_seqAndRt)
    {
        TransformationDescription::DataPoints trafo_data_tmp;
        auto trafoit = trafo_tmp.back().getDataPoints().begin();
        auto mapit = mapsit.begin();
        while (trafoit != trafo_tmp.back().getDataPoints().end() && mapit != mapsit.end())
        {
            if (trafoit->note < mapit->first)
            {
                ++trafoit;
            }
            else if (trafoit->note > mapit->first)
            {
                ++mapit;
            }
            else
            {
                // TODO: check problems with outliers
                for (auto & rt : mapit->second)
                {
                    TransformationDescription::DataPoint point(rt, trafoit->second, trafoit->note);
                    trafo_data_tmp.push_back(point);
                }
                ++trafoit;
                ++mapit;
            }
        }
        TransformationDescription map_trafo = TransformationDescription(trafo_data_tmp);
        map_trafo.fitModel(model_type, model_params);
        transformations.push_back(map_trafo);
    }
}

int main(int argc, const char** argv)
{
  TOPPMapAlignerTree tool;
  return tool.main(argc, argv);
}

/// @endcond
