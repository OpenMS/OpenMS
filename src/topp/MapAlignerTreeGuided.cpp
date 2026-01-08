// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julia Thueringer $
// $Authors: Julia Thueringer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h> // to print newick tree on cml


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MapAlignerTreeGuided MapAlignerTreeGuided

@brief Corrects retention time distortions between maps, using information from peptides identified in different maps.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN= "middle" ROWSPAN=2> &rarr; MapAlignerTreeGuided &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper @n (or any other source of FDR-filtered
featureXMLs) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeledKD or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
        </tr>
    </table>
</CENTER>


This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and
distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g.
before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with
peptide identifications obtained in a different run.

All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to
this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations
right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the
transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to
arbitrary files with the @ref TOPP_MapRTTransformer tool.

The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and
consequently what types of data they can be applied to. The alignment algorithm implemented here is based on peptide
identifications and applicable to annotated featureXML files. It finds peptide sequences that each pair of input files
have in common, uses them as points of correspondence between the inputs and to evaluate the distances between the
maps for hierarchical clustering. Tree based, the alignment of each cluster pair is performed with the method align() of the
@ref OpenMS::MapAlignmentAlgorithmIdentification. For more details and algorithm-specific parameters (set in the INI file)
see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmTreeGuided "algorithm documentation".

@see @ref TOPP_MapAlignerIdentification @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapRTTransformer

Note that alignment is based on the sequence including modifications, thus an exact match is required. I.e., a peptide with oxidised methionine will not be matched to its unmodified version. This behavior is generally desired since (some) modifications can cause retention time shifts.

Also note that convex hulls are removed for alignment and are therefore missing in the output files.

Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested with the "b_spline" model. The different available models are:
- @ref OpenMS::TransformationModelLinear "linear": Linear model.
- @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
- @ref OpenMS::TransformationModelLowess "lowess": Local regression (non-linear).
- @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

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
          TOPPMapAlignerBase("MapAlignerTreeGuided", "Tree guided correction of retention time distortions between maps.")
  {
  }

private:
  template <typename MapType>
  void loadInputMaps_(vector<MapType>& maps, StringList& ins, FileHandler& fxml_file)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
    progresslogger.startProgress(0, ins.size(), "loading input files");
    for (Size i = 0; i < ins.size(); ++i)
    {
      progresslogger.setProgress(i);
      fxml_file.loadFeatures(ins[i], maps[i], {FileTypes::FEATUREXML});
    }
    progresslogger.endProgress();
  }

  void storeFeatureXMLs_(vector<FeatureMap>& feature_maps, const StringList& out_files, FileHandler& fxml_file)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
    progresslogger.startProgress(0, feature_maps.size(), "writing output files");
    for (Size i = 0; i < out_files.size(); ++i)
    {
      progresslogger.setProgress(i);
      // annotate output with data processing info
      addDataProcessing_(feature_maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));
      fxml_file.storeFeatures(out_files[i], feature_maps[i], {FileTypes::FEATUREXML});
    }
    progresslogger.endProgress();
  }

  void storeTransformationDescriptions_(const vector<TransformationDescription>& transformations, StringList& trafos)
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
    progresslogger.startProgress(0, trafos.size(),
                                 "writing transformation files");
    for (Size i = 0; i < trafos.size(); ++i)
    {
      progresslogger.setProgress(i);
      FileHandler().storeTransformations(trafos[i], transformations[i], {FileTypes::TRANSFORMATIONXML});
    }
    progresslogger.endProgress();
  }

  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlagsMapAligners_("featureXML",
                                                            REF_NONE);
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerStringOption_("copy_data", "String", "true", "Copy data (faster, more memory required) or reload data (slower, less memory required) when aligning many files.", false, false);
    setValidStrings_("copy_data", {"true","false"});
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmTreeGuided algo;
      return algo.getParameters();
    }
    return Param(); // shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes ret = checkParameters_();
    if (ret != EXECUTION_OK)
    {
      return ret;
    }
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    StringList in_files = getStringList_("in");
    StringList out_files = getStringList_("out");
    StringList out_trafos = getStringList_("trafo_out");
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    Size in_files_size = in_files.size();
    FileHandler fxml_file;
    // define here because needed to load and store
    FeatureFileOptions param = fxml_file.getFeatOptions();
    // to save memory don't load convex hulls and subordinates
    param.setLoadSubordinates(false);
    param.setLoadConvexHull(false);
    fxml_file.setFeatOptions(param);

    vector<FeatureMap> feature_maps(in_files_size);
    loadInputMaps_(feature_maps, in_files, fxml_file);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // constructing tree
    vector<vector<double>> maps_ranges(in_files_size);  // to save ranges for alignment (larger rt_range -> reference)
    std::vector<BinaryTreeNode> tree;    // to construct tree with pearson coefficient
    MapAlignmentAlgorithmTreeGuided algoTree;
    Param algo_params = getParam_().copy("algorithm:", true);
    algoTree.setParameters(algo_params);
    OpenMS::MapAlignmentAlgorithmTreeGuided::buildTree(feature_maps, tree, maps_ranges);

    // print tree
    ClusterAnalyzer ca;
    OPENMS_LOG_INFO << "  Alignment follows Newick tree: " << ca.newickTree(tree, true) << endl;

    // alignment
    vector<Size> trafo_order;
    FeatureMap map_transformed;
    // depending on the selected parameter, the input data for the alignment are copied or reloaded after alignment
    if (getStringOption_("copy_data") == "true")
    {
      vector<FeatureMap> copied_maps = feature_maps;
      algoTree.treeGuidedAlignment(tree, copied_maps, maps_ranges, map_transformed, trafo_order);
    }
    else
    {
      algoTree.treeGuidedAlignment(tree, feature_maps, maps_ranges, map_transformed, trafo_order);
      // load() of FeatureXMLFile clears featureMap, so we don't have to care
      loadInputMaps_(feature_maps, in_files, fxml_file);
    }

    //-------------------------------------------------------------
    // generating output
    //-------------------------------------------------------------
    vector<TransformationDescription> transformations(in_files_size); // for trafo_out
    algoTree.computeTrafosByOriginalRT(feature_maps, map_transformed, transformations, trafo_order);
    OpenMS::MapAlignmentAlgorithmTreeGuided::computeTransformedFeatureMaps(feature_maps, transformations);

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

int main(int argc, const char** argv)
{
  TOPPMapAlignerTreeGuided tool;
  return tool.main(argc, argv);
}

/// @endcond

