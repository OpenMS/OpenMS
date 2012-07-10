// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_MapAlignerIdentification MapAlignerIdentification

        @brief Corrects retention time distortions between maps, using information from peptides identified in different maps.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ MapAlignerIdentification \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter @n (or another search engine adapter) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMerger </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
        </tr>
    </table>
</CENTER>

    This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

    All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

    The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types of data they can be applied to. The alignment algorithm implemented here is based on peptide identifications, and thus applicable to files containing peptide IDs (idXML, annotated featureXML/consensusXML). It finds peptide sequences that different input files have in common and uses them as points of correspondence between the inputs. For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmIdentification "algorithm documentation".

    @see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum @ref TOPP_MapRTTransformer

    Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested mostly with the "b_spline" model. The different available models are:
    - @ref OpenMS::TransformationModelLinear "linear": Linear model.
    - @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
    - @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

    The following parameters control the modeling of RT transformations (they can be set in the "model" section of the INI file):
    @htmlinclude OpenMS_MapAlignerIdentificationModel.parameters @n

    <B>The command line parameters of this tool are:</B> @n
    @verbinclude TOPP_MapAlignerIdentification.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_MapAlignerIdentification.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerIdentification :
  public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerIdentification() :
    TOPPMapAlignerBase("MapAlignerIdentification", "Corrects retention time distortions between maps based on common peptide identifications.")
  {
  }

private:
  template <typename TMapType, typename TFileType>
  void loadInitialMaps(std::vector<TMapType> & maps, StringList & ins, TFileType & inputPutFile)
  {
    // custom progresslogger for this task
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, ins.size(), "loading input files");

    for (Size i = 0; i < ins.size(); ++i)
    {
      progresslogger.setProgress(i);
      inputPutFile.load(ins[i], maps[i]);
    }

    progresslogger.endProgress();
  }

  /// helper function to avoid code duplication between consensus and feautreXML storage operations
  template <typename TMapType, typename TFileType>
  void storeTransformedMaps(std::vector<TMapType> & maps, StringList & outs, TFileType & outPutFile)
  {
    // custom progresslogger for this task
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, outs.size(), "writing output files");

    for (Size i = 0; i < outs.size(); ++i)
    {
      progresslogger.setProgress(i);

      //annotate output with data processing info
      addDataProcessing_(maps[i], getProcessingInfo_(DataProcessing::ALIGNMENT));

      outPutFile.store(outs[i], maps[i]);
    }
    progresslogger.endProgress();
  }

  void saveTransformationDescriotions(const std::vector<TransformationDescription> & transformations,
                                      StringList & trafos)
  {
    // custom progresslogger for this task
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, trafos.size(), "writing transformation files");

    for (Size i = 0; i < transformations.size(); ++i)
    {
      TransformationXMLFile().store(trafos[i], transformations[i]);
    }

    progresslogger.endProgress();
  }

  void registerOptionsAndFlags_()
  {
    String formats = "featureXML,consensusXML,idXML";
    TOPPMapAlignerBase::registerOptionsAndFlags_(formats, true);
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  Param getSubsectionDefaults_(const String & section) const
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return getModelDefaults("b_spline");
    }

    // shouldn't happen
    return Param();
  }

  ExitCodes main_(int, const char **)
  {
    MapAlignmentAlgorithmIdentification algorithm;
    handleReference_(&algorithm);

    ExitCodes returnCode = initialize_(&algorithm);
    if (returnCode != EXECUTION_OK) return returnCode;

    // handle in- and output files
    StringList inputFiles = getStringList_("in");
    StringList outputFiles = getStringList_("out");
    StringList trafoFiles = getStringList_("trafo_out");
    FileTypes::Type in_type = FileHandler::getType(inputFiles[0]);

    // find model parameters
    Param modelParams = getParam_().copy("model:", true);
    String modelType = modelParams.getValue("type");
    modelParams = modelParams.copy(modelType + ":", true);

    // create transformations vector
    std::vector<TransformationDescription> transformations;

    if (in_type == FileTypes::FEATUREXML)
    {
      // load input
      std::vector<FeatureMap<> > featureMaps(inputFiles.size());
      FeatureXMLFile fxmlFile;

      // no need to store featureXML, thus we can load only minimum required information
      if (outputFiles.size() == 0)
      {
        fxmlFile.getOptions().setLoadConvexHull(false);
        fxmlFile.getOptions().setLoadSubordinates(false);
      }

      // load maps
      loadInitialMaps(featureMaps, inputFiles, fxmlFile);

      algorithm.alignFeatureMaps(featureMaps, transformations);

      if (modelType != "none")
      {
        algorithm.fitModel(modelType, modelParams, transformations);
      }

      MapAlignmentTransformer::transformFeatureMaps(featureMaps, transformations);

      storeTransformedMaps(featureMaps, outputFiles, fxmlFile);
    }
    //-------------------------------------------------------------
    // perform consensus alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::CONSENSUSXML)
    {
      // load input
      std::vector<ConsensusMap> cons_maps(inputFiles.size());
      ConsensusXMLFile f;

      // load maps
      loadInitialMaps(cons_maps, inputFiles, f);

      algorithm.alignConsensusMaps(cons_maps, transformations);

      if (modelType != "none")
      {
        algorithm.fitModel(modelType, modelParams, transformations);
      }
      MapAlignmentTransformer::transformConsensusMaps(cons_maps, transformations);

      storeTransformedMaps(cons_maps, outputFiles, f);
    }
    //-------------------------------------------------------------
    // perform peptide alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::IDXML)
    {
      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);

      // load input
      std::vector<std::vector<ProteinIdentification> > protein_ids_vec(inputFiles.size());
      std::vector<std::vector<PeptideIdentification> > peptide_ids_vec(inputFiles.size());

      IdXMLFile f;

      progresslogger.startProgress(0, inputFiles.size(), "loading input files");
      for (Size i = 0; i < inputFiles.size(); ++i)
      {
        progresslogger.setProgress(i);
        f.load(inputFiles[i], protein_ids_vec[i], peptide_ids_vec[i]);
      }
      progresslogger.endProgress();

      algorithm.alignPeptideIdentifications(peptide_ids_vec, transformations);

      if (modelType != "none")
      {
        algorithm.fitModel(modelType, modelParams, transformations);
      }

      MapAlignmentTransformer::transformPeptideIdentifications(peptide_ids_vec,
                                                               transformations);

      // write output
      progresslogger.startProgress(0, outputFiles.size(), "writing output files");
      for (Size i = 0; i < outputFiles.size(); ++i)
      {
        progresslogger.setProgress(i);
        f.store(outputFiles[i], protein_ids_vec[i], peptide_ids_vec[i]);
      }
      progresslogger.endProgress();
    }
    else
    {
      // TODO can this really happen? I think it is tested above. Otherwise
      // throw an appropriate exception?
      return ILLEGAL_PARAMETERS;
    }

    if (!trafoFiles.empty())
    {
      saveTransformationDescriotions(transformations, trafoFiles);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPMapAlignerIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
