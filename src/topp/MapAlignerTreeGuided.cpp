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
// $Authors: $
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

  @brief Corrects retention time distortions between maps by following a tree constructed with SingleLinkage

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ MapAlignerTreeGuided \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_XTandemAdapter @n (or another search engine adapter) </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinkerUnlabeledKD or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
        </tr>
    </table>
</CENTER>

    This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

    All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

    The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types of data they can be applied to. The alignment algorithm implemented here is based on peptide identifications and applicable to annotated featureXML files. It finds peptide sequences that different input files have in common and uses them as points of correspondence between the inputs. For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmIdentification "algorithm documentation".

    @see @ref TOPP_MapAlignerIdentification @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum @ref TOPP_MapRTTransformer

		Note that alignment is based on the sequence including modifications, thus an exact match is required. I.e., a peptide with oxidised methionine will not be matched to its unmodified version. This behavior is generally desired since (some) modifications can cause retention time shifts.

    Also note that convex hulls are removed for alignment and are therefore missing in the output files.

    Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested with the "b_spline" model. The different available models are:
    - @ref OpenMS::TransformationModelLinear "linear": Linear model.
    - @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
    - @ref OpenMS::TransformationModelLowess "lowess": Local regression (non-linear).
    - @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

    The following parameters control the modeling of RT transformations (they can be set in the "model" section of the INI file):
    @htmlinclude OpenMS_MapAlignerTreeGuidedModel.parameters @n

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B> @n
    @verbinclude TOPP_MapAlignerTreeGuided.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MapAlignerTreeGuided.html
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

      return ((1 - pearson_val) * intercept_size / union_size);
    }
  }; // end of PeptideIdentificationsPearsonDifference

  // function declarations
  template <typename MapType>
  void loadInputMaps_(vector<MapType>& maps, StringList& ins, FeatureXMLFile& fxml_file);
  static void getPeptideSequences_(const vector<PeptideIdentification>& peptides, SeqAndRTList& peptide_rts, vector<double>& rts_tmp);
  static void extractSeqAndRt_(const vector<FeatureMap>& feature_maps, vector<SeqAndRTList>& maps_seq_and_rt, vector<vector<double>>& maps_ranges);
  static void buildTree_(vector<FeatureMap>& feature_maps, std::vector<BinaryTreeNode>& tree, vector<vector<double>>& maps_ranges);
  void treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree, vector<FeatureMap> feature_maps_transformed,
          vector<vector<double>> &maps_ranges, FeatureMap& map_transformed, vector<Size>& trafo_order, const Param& model_params, const String& model_type);
  static void computeTransformationsByOrigInFeature_(vector<FeatureMap>& feature_maps, FeatureMap& map_transformed,
          vector<TransformationDescription>& transformations, const vector<Size>& trafo_order, const Param& model_params,
          const String& model_type);
  void computeTransformedFeatureMaps_(vector<FeatureMap>& feature_maps, const vector<TransformationDescription>& transformations);
  void storeFeatureXMLs_(const vector<FeatureMap>& feature_maps, const StringList& out_files, FeatureXMLFile& fxml_file);
  void storeTransformationDescriptions_(const vector<TransformationDescription>& transformations, StringList& trafos);


  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlags_("featureXML", REF_FLEXIBLE);
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

    vector<vector<double>> maps_ranges(in_files_size);  // to save ranges for alignment (larger rt_range -> reference)
    std::vector<BinaryTreeNode> tree;    // to construct tree with pearson coefficient
    buildTree_(feature_maps, tree, maps_ranges);

    // print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "  Alignment follows tree: " << ca.newickTree(tree) << endl;

    vector<Size> trafo_order;
    FeatureMap map_transformed;
    // alignment uses copy of feature_maps, returning result within map_transformed (to keep original data)
    treeGuidedAlignment_(tree, feature_maps, maps_ranges, map_transformed, trafo_order, model_params, model_type);

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

void TOPPMapAlignerTreeGuided::extractSeqAndRt_(const vector<FeatureMap>& feature_maps,
                                                vector<SeqAndRTList>& maps_seq_and_rt, vector<vector<double>>& maps_ranges)
{
  for (Size i = 0; i < feature_maps.size(); ++i)
  {
    vector<double> rts_tmp(feature_maps[i].size());
    for (auto feature_it = feature_maps[i].begin(); feature_maps[i].end() != feature_it; ++feature_it)
    {
      if (!feature_it->getPeptideIdentifications().empty()) {
        getPeptideSequences_(feature_it->getPeptideIdentifications(), maps_seq_and_rt[i], rts_tmp);
      }
    }
    sort(rts_tmp.begin(), rts_tmp.end());

    maps_ranges[i] = rts_tmp;
    rts_tmp.clear();
  }
}

void TOPPMapAlignerTreeGuided::buildTree_(vector<FeatureMap>& feature_maps, std::vector<BinaryTreeNode>& tree, vector<vector<double>>& maps_ranges)
{
  vector<SeqAndRTList> maps_seq_and_rt(feature_maps.size());
  extractSeqAndRt_(feature_maps, maps_seq_and_rt, maps_ranges);
  PeptideIdentificationsPearsonDistance pep_dist;
  SingleLinkage sl;
  DistanceMatrix<float> dist_matrix;
  ClusterHierarchical ch;
  ch.cluster<SeqAndRTList, PeptideIdentificationsPearsonDistance>(maps_seq_and_rt, pep_dist, sl, tree, dist_matrix);
}

void TOPPMapAlignerTreeGuided::treeGuidedAlignment_(const std::vector<BinaryTreeNode> &tree, vector<FeatureMap> feature_maps_transformed,
        vector<vector<double>> &maps_ranges, FeatureMap& map_transformed,
        vector<Size>& trafo_order, const Param& model_params, const String& model_type)
{
  Size last_trafo = 0;  // to get final transformation order from map_sets
  vector<TransformationDescription> transformations_align;  // temporary for aligner output

  // helper to memorize rt transformation order
  vector<vector<Size>> map_sets(feature_maps_transformed.size());
  for (Size i = 0; i < feature_maps_transformed.size(); ++i)
  {
    map_sets[i].push_back(i);
  }

  MapAlignmentAlgorithmIdentification algorithm;
  Param algo_params = getParam_().copy("align_algorithm:", true);
  algorithm.setParameters(algo_params);
  algorithm.setLogType(log_type_);

  bool use_ranges = getStringOption_("use_ranges").toQString() == "true";
  Size ref;
  Size to_transform;

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
      double left_range = maps_ranges[node.left_child][maps_ranges[node.left_child].size()*0.9] - maps_ranges[node.left_child][maps_ranges[node.left_child].size()*0.1];
      double right_range = maps_ranges[node.right_child][maps_ranges[node.right_child].size()*0.9] - maps_ranges[node.right_child][maps_ranges[node.right_child].size()*0.1];

      if (left_range > right_range) {
        ref = node.left_child;
        to_transform = node.right_child;
      } else {
        ref = node.right_child;
        to_transform = node.left_child;
      }
      std::vector<double> tmp;
      std::merge(maps_ranges[node.right_child].begin(), maps_ranges[node.right_child].end(), maps_ranges[node.left_child].begin(), maps_ranges[node.left_child].end(), std::back_inserter(tmp));
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
    OPENMS_LOG_INFO << "use internal reference for alignment"  << endl;
    algorithm.align(to_align, transformations_align, 1);
    // transform retention times of non-identity for next iteration
    transformations_align[0].fitModel(model_type, model_params);
    MapAlignmentTransformer::transformRetentionTimes(feature_maps_transformed[to_transform],
                                                     transformations_align[0], true);

    // combine aligned maps, store in both, because tree always calls smaller number
    // also possible: feature_maps_transformed[smallerNumber] = ..[ref]+..[to_transform]
    // or use pointer
    feature_maps_transformed[ref] += feature_maps_transformed[to_transform];
    feature_maps_transformed[ref].updateRanges();
    feature_maps_transformed[to_transform] = feature_maps_transformed[ref];

    // update order of alignment for both aligned maps
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

