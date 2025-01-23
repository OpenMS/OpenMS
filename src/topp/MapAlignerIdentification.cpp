// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/OMSFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MapAlignerIdentification MapAlignerIdentification

@brief Corrects retention time distortions between maps, using information from peptides identified in different maps.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> &rarr; MapAlignerIdentification &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
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

Reference:\n
Weisser <em>et al.</em>: <a href="https://doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types of data they can be applied to. The alignment algorithm implemented here is based on peptide identifications, and thus applicable to files containing peptide IDs (idXML, annotated featureXML/consensusXML). It finds peptide sequences that different input files have in common and uses them as points of correspondence between the inputs. For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmIdentification "algorithm documentation".

@see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapRTTransformer

Note that alignment is based on the sequence including modifications, thus an exact match is required. I.e., a peptide with oxidised methionine will not be matched to its unmodified version. This behavior is generally desired since (some) modifications can cause retention time shifts.

Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested mostly with the "b_spline" model. The different available models are:
- @ref OpenMS::TransformationModelLinear "linear": Linear model.
- @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
- @ref OpenMS::TransformationModelLowess "lowess": Local regression (non-linear).
- @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

The following parameters control the modeling of RT transformations (they can be set in the "model" section of the INI file):
@htmlinclude OpenMS_MapAlignerIdentificationModel.parameters @n

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

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
  template <typename MapType, typename FileType>
  void loadInitialMaps_(vector<MapType>& maps, StringList& ins,
                        FileType& input_file)
  {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
    progresslogger.startProgress(0, ins.size(), "loading input files");
    for (Size i = 0; i < ins.size(); ++i)
    {
      progresslogger.setProgress(i);
      input_file.load(ins[i], maps[i]);
    }
    progresslogger.endProgress();
  }

  // helper function to avoid code duplication between consensusXML and
  // featureXML storage operations:
  template <typename MapType, typename FileType>
  void storeTransformedMaps_(vector<MapType>& maps, StringList& outs,
                             FileType& output_file)
  {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, outs.size(), "writing output files");
    for (Size i = 0; i < outs.size(); ++i)
    {
      progresslogger.setProgress(i);
      // annotate output with data processing info:
      addDataProcessing_(maps[i],
                         getProcessingInfo_(DataProcessing::ALIGNMENT));
      output_file.store(outs[i], maps[i]);
    }
    progresslogger.endProgress();
  }

  template <typename DataType>
  void performAlignment_(MapAlignmentAlgorithmIdentification& algorithm,
                         vector<DataType>& data,
                         vector<TransformationDescription>& transformations,
                         Int reference_index)
  {
    // find model parameters:
    Param model_params = getParam_().copy("model:", true);
    String model_type = model_params.getValue("type").toString();

    try
    {
      algorithm.align(data, transformations, reference_index);
    }
    catch (Exception::MissingInformation& err)
    {
      if (getFlag_("force"))
      {
        OPENMS_LOG_ERROR
          << "Error: alignment failed. Details:\n" << err.what()
          << "\nSince 'force' is set, processing will continue using 'identity' transformations."
          << endl;
        model_type = "identity";
        transformations.resize(data.size());
      }
      else throw;
    }

    if (model_type != "none")
    {
      model_params = model_params.copy(model_type + ":", true);
      for (TransformationDescription& tra : transformations)
      {
        tra.fitModel(model_type, model_params);
      }
    }
  }

  template <typename DataType>
  void applyTransformations_(vector<DataType>& data,
    const vector<TransformationDescription>& transformations)
  {
    bool store_original_rt = getFlag_("store_original_rt");
    for (Size i = 0; i < data.size(); ++i)
    {
      MapAlignmentTransformer::transformRetentionTimes(
        data[i], transformations[i], store_original_rt);
    }
  }

  void storeTransformationDescriptions_(const vector<TransformationDescription>&
                                        transformations, StringList& trafos)
  {
    OPENMS_PRECONDITION(transformations.size() == trafos.size(), "Transformation descriptions and list of transformation files need to be equal.");
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, trafos.size(),
                                 "writing transformation files");
    OPENMS_LOG_INFO << "Writing " << transformations.size() << " transformations " << " to " << trafos.size() << " files.";
    for (Size i = 0; i < transformations.size(); ++i)
    {
      FileHandler().storeTransformations(trafos[i], transformations[i], {FileTypes::TRANSFORMATIONXML});
    }
    progresslogger.endProgress();
  }

  Int getReference_(MapAlignmentAlgorithmIdentification& algorithm)
  {
    // consistency of reference parameters has already been checked via
    // "TOPPMapAlignerBase::checkParameters_"

    Size reference_index = getIntOption_("reference:index");
    String reference_file = getStringOption_("reference:file");

    if (!reference_file.empty())
    {
      FileTypes::Type filetype = FileHandler::getType(reference_file);
      switch (filetype)
      {
      case FileTypes::MZML:
      {
        PeakMap experiment;
        FileHandler().loadExperiment(reference_file, experiment, {FileTypes::MZML});
        algorithm.setReference(experiment);
      }
      break;
      case FileTypes::FEATUREXML:
      {
        FeatureMap features;
        FileHandler().loadFeatures(reference_file, features);
        algorithm.setReference(features);
      }
      break;
      case FileTypes::CONSENSUSXML:
      {
        ConsensusMap consensus;
        FileHandler().loadConsensusFeatures(reference_file, consensus);
        algorithm.setReference(consensus);
      }
      break;
      case FileTypes::IDXML:
      {
        vector<ProteinIdentification> proteins;
        vector<PeptideIdentification> peptides;
        FileHandler().loadIdentifications(reference_file, proteins, peptides);
        algorithm.setReference(peptides);
      }
      break;
      case FileTypes::OMS:
      {
        IdentificationData id_data;
        OMSFile().load(reference_file, id_data);
        algorithm.setReference(id_data);
      }
      break;
      default: // to avoid compiler warnings
        throw Exception::WrongParameterType(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION,
                                            "reference:file");
      }
    }

    return Int(reference_index) - 1; // internally, we count from zero
  }

  void registerOptionsAndFlags_() override
  {
    String formats = "featureXML,consensusXML,idXML,oms";
    TOPPMapAlignerBase::registerOptionsAndFlagsMapAligners_(formats, REF_FLEXIBLE);
    // TODO: potentially move to base class so every aligner has to support design
    registerInputFile_("design", "<file>", "", "Input file containing the experimental design", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    registerFlag_("store_original_rt", "Store the original retention times (before transformation) as meta data in the output?");

    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return MapAlignerBase::getModelDefaults("b_spline");
    }

    return Param(); // this shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes return_code = TOPPMapAlignerBase::checkParameters_();
    if (return_code != EXECUTION_OK) return return_code;

    // set up alignment algorithm:
    MapAlignmentAlgorithmIdentification algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    Int reference_index = getReference_(algorithm);

    // handle in- and output files:
    StringList input_files = getStringList_("in");
    if (input_files.size() == 1)
    {
      OPENMS_LOG_WARN << "Only one file provided as input to MapAlignerIdentification." << std::endl;
    }   

    StringList output_files = getStringList_("out");
    StringList trafo_files = getStringList_("trafo_out");
    FileTypes::Type in_type = FileHandler::getType(input_files[0]);

    vector<TransformationDescription> transformations;

    switch (in_type)
    {
    //-------------------------------------------------------------
    // perform feature alignment
    //-------------------------------------------------------------
    case FileTypes::FEATUREXML:
    {
      vector<FeatureMap> feature_maps(input_files.size());
      FeatureXMLFile fxml_file;
      if (output_files.empty())
      {
        // store only transformation descriptions, not transformed data =>
        // we can load only minimum required information:
        fxml_file.getOptions().setLoadConvexHull(false);
        fxml_file.getOptions().setLoadSubordinates(false);
      }
      loadInitialMaps_(feature_maps, input_files, fxml_file);

      //-------------------------------------------------------------
      // extract (optional) fraction identifiers and associate with featureXMLs
      //-------------------------------------------------------------
      String design_file = getStringOption_("design");

      // determine map of fractions to runs
      map<unsigned, vector<String>> frac2files;

      // TODO: check if can be put in common helper function
      if (!design_file.empty())
      {
        // parse design file and determine fractions
        ExperimentalDesign ed = ExperimentalDesignFile::load(design_file,
                                                             false);

        // determine if design defines more than one fraction (note: fraction and run IDs are one-based)
        frac2files = ed.getFractionToMSFilesMapping();

        // check if all fractions have the same number of MS runs associated
        if (!ed.sameNrOfMSFilesPerFraction())
        {
          writeLogError_("Error: Number of runs must match for every fraction!");
          return ILLEGAL_PARAMETERS;
        }
      }
      else // no design file given
      {
        for (Size i = 0; i != input_files.size(); ++i)
        {
          // TODO: read proper MS file name from meta data
          frac2files[1].push_back("file" + String(i)); // associate each file with fraction 1
        }
      }

      // TODO: check and handle if featureXML order differs from run order

      // perform fraction-based alignment
      if (frac2files.size() == 1) // group one fraction
      {
        performAlignment_(algorithm, feature_maps, transformations,
                          reference_index);
        applyTransformations_(feature_maps, transformations);
      }
      else // group multiple fractions
      {
        for (Size i = 1; i <= frac2files.size(); ++i)
        {
          vector<FeatureMap> fraction_maps;
          vector<TransformationDescription> fraction_transformations;

          size_t n_fractions = frac2files.size();

          // TODO FRACTIONS: determine map index based on annotated MS files (getPrimaryMSRuns())
          for (size_t feature_map_index = 0; feature_map_index != n_fractions;
               ++feature_map_index)
          {
            fraction_maps.push_back(feature_maps[feature_map_index]);
          }
          performAlignment_(algorithm, fraction_maps, fraction_transformations,
                            reference_index);
          applyTransformations_(fraction_maps, fraction_transformations);

          // copy into transformations and feature maps
          transformations.insert(transformations.end(),
                                 fraction_transformations.begin(),
                                 fraction_transformations.end());

          Size f = 0;
          for (size_t feature_map_index = 0; feature_map_index != n_fractions;
               ++feature_map_index, ++f)
          {
            feature_maps[feature_map_index].swap(fraction_maps[f]);
          }
        }
      }

      if (!output_files.empty())
      {
        storeTransformedMaps_(feature_maps, output_files, fxml_file);
      }
    }
    break;

    //-------------------------------------------------------------
    // perform consensus alignment
    //-------------------------------------------------------------
    case FileTypes::CONSENSUSXML:
    {
      std::vector<ConsensusMap> consensus_maps(input_files.size());
      ConsensusXMLFile cxml_file;
      loadInitialMaps_(consensus_maps, input_files, cxml_file);

      performAlignment_(algorithm, consensus_maps, transformations,
                        reference_index);
      applyTransformations_(consensus_maps, transformations);

      if (!output_files.empty())
      {
        storeTransformedMaps_(consensus_maps, output_files, cxml_file);
      }
    }
    break;

    //-------------------------------------------------------------
    // perform peptide alignment
    //-------------------------------------------------------------
    case FileTypes::IDXML:
    {
      vector<vector<ProteinIdentification>> protein_ids(input_files.size());
      vector<vector<PeptideIdentification>> peptide_ids(input_files.size());
      FileHandler idxml_file;
      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);
      progresslogger.startProgress(0, input_files.size(),
                                   "loading input files");
      for (Size i = 0; i < input_files.size(); ++i)
      {
        progresslogger.setProgress(i);
        idxml_file.loadIdentifications(input_files[i], protein_ids[i], peptide_ids[i], {FileTypes::IDXML});
      }
      progresslogger.endProgress();

      performAlignment_(algorithm, peptide_ids, transformations,
                        reference_index);
      applyTransformations_(peptide_ids, transformations);

      if (!output_files.empty())
      {
        progresslogger.startProgress(0, output_files.size(),
                                     "writing output files");
        for (Size i = 0; i < output_files.size(); ++i)
        {
          progresslogger.setProgress(i);
          idxml_file.storeIdentifications(output_files[i], protein_ids[i], peptide_ids[i], {FileTypes::IDXML});
        }
        progresslogger.endProgress();
      }
    }
    break;

    //-------------------------------------------------------------
    // perform spectrum match alignment
    //-------------------------------------------------------------
    case FileTypes::OMS:
    {
      vector<IdentificationData> id_data(input_files.size());
      OMSFile oms_file;
      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);
      progresslogger.startProgress(0, input_files.size(),
                                   "loading input files");
      for (Size i = 0; i < input_files.size(); ++i)
      {
        progresslogger.setProgress(i);
        oms_file.load(input_files[i], id_data[i]);
      }
      progresslogger.endProgress();

      // add data processing information:
      DateTime processing_time = DateTime::now(); // use same for each file
      IdentificationData::ProcessingSoftware sw(toolName_(), version_);
      if (test_mode_) sw.setVersion("test");
      String reference_file = getStringOption_("reference:file");
      for (IdentificationData& id : id_data)
      {
        IdentificationData::ProcessingSoftwareRef sw_ref =
          id.registerProcessingSoftware(sw);
        IdentificationData::ProcessingStep step(sw_ref);
        for (const String& input_file : input_files)
        {
          IdentificationData::InputFileRef ref =
            id.registerInputFile(IdentificationData::InputFile(input_file));
          step.input_file_refs.push_back(ref);
        }
        if (!reference_file.empty())
        {
          IdentificationData::InputFileRef ref =
            id.registerInputFile(IdentificationData::InputFile(reference_file));
          step.input_file_refs.push_back(ref);
        }
        step.date_time = processing_time;
        step.actions.insert(DataProcessing::ALIGNMENT);
        id.registerProcessingStep(step);
      }

      performAlignment_(algorithm, id_data, transformations, reference_index);
      applyTransformations_(id_data, transformations);

      if (!output_files.empty())
      {
        progresslogger.startProgress(0, output_files.size(),
                                     "writing output files");
        for (Size i = 0; i < output_files.size(); ++i)
        {
          progresslogger.setProgress(i);
          oms_file.store(output_files[i], id_data[i]);
        }
        progresslogger.endProgress();
      }
    }
    break;

    default: // to avoid compiler warnings
      throw Exception::WrongParameterType(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, "in");
    }

    if (!trafo_files.empty())
    {
      storeTransformationDescriptions_(transformations, trafo_files);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPMapAlignerIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
