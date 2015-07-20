// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_MAPALIGNERBASE_H
#define OPENMS_APPLICATIONS_MAPALIGNERBASE_H

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_MapAlignerBase MapAlignerBase

    @brief Base class for different MapAligner TOPP tools.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

class TOPPMapAlignerBase :
  public TOPPBase
{

public:
  TOPPMapAlignerBase(String name, String description, bool official = true) :
    TOPPBase(name, description, official)
  {
  }

  // "public" so it can be used in DefaultParamHandlerDocumenter to get docu
  static Param getModelDefaults(const String& default_model)
  {
    Param params;
    params.setValue("type", default_model, "Type of model");
    // TODO: avoid referring to each TransformationModel subclass explicitly
    StringList model_types = ListUtils::create<String>("linear,b_spline,interpolated");
    if (!ListUtils::contains(model_types, default_model))
    {
      model_types.insert(model_types.begin(), default_model);
    }
    params.setValidStrings("type", model_types);

    Param model_params;
    TransformationModelLinear::getDefaultParameters(model_params);
    params.insert("linear:", model_params);
    params.setSectionDescription("linear", "Parameters for 'linear' model");
    TransformationModelBSpline::getDefaultParameters(model_params);
    params.insert("b_spline:", model_params);
    params.setSectionDescription("b_spline", "Parameters for 'b_spline' model");
    TransformationModelInterpolated::getDefaultParameters(model_params);
    params.insert("interpolated:", model_params);
    params.setSectionDescription("interpolated",
                                 "Parameters for 'interpolated' model");
    return params;
  }

protected:
  void registerOptionsAndFlags_(const String& file_formats, const bool add_reference = false)
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blanks (all must have the same file type)", true);
    setValidFormats_("in", ListUtils::create<String>(file_formats));
    registerOutputFileList_("out", "<files>", StringList(), "Output files separated by blanks. Either 'out' or 'trafo_out' has to be provided. They can be used together.", false);
    setValidFormats_("out", ListUtils::create<String>(file_formats));
    registerOutputFileList_("trafo_out", "<files>", StringList(), "Transformation output files separated by blanks. Either 'out' or 'trafo_out' has to be provided. They can be used together.", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    addEmptyLine_();
    if (add_reference)
    {
      registerTOPPSubsection_("reference", "Options to define a reference file (use either 'file' or 'index', not both; if neither is given 'index' is used).");
      registerInputFile_("reference:file", "<file>", "", "File to use as reference (same file format as input files required)", false);
      setValidFormats_("reference:file", ListUtils::create<String>(file_formats));
      registerIntOption_("reference:index", "<number>", 0, "Use one of the input files as reference ('1' for the first file, etc.).\nIf '0', no explicit reference is set - the algorithm will select a reference.", false);
      setMinInt_("reference:index", 0);
    }
  }

  /// deprecated? (not used in PoseClustering... and moved to initialize_() )
  void handleReference_(MapAlignmentAlgorithm* alignment)
  {
    // note: this function is in the base class to avoid code duplication, but
    // it only makes sense for some derived classes - don't call the function
    // in a class that doesn't support a reference!

    // check reference parameters:
    Size reference_index = getIntOption_("reference:index");
    String reference_file = getStringOption_("reference:file");
    if (reference_index > getStringList_("in").size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'reference:index' must not be higher than the number of input files");
    }
    if (reference_index && !reference_file.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'reference:index' and 'reference:file' cannot be used together");
    }

    // pass the reference parameters on to the algorithm:
    alignment->setReference(reference_index, reference_file);
  }

  ExitCodes initialize_(MapAlignmentAlgorithm* alignment, bool check_ref = false)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    StringList ins = getStringList_("in");
    StringList outs = getStringList_("out");
    StringList trafos = getStringList_("trafo_out");

    //-------------------------------------------------------------
    // check for valid input
    //-------------------------------------------------------------
    // check whether some kind of output file is given:
    if (outs.empty() && trafos.empty())
    {
      writeLog_("Error: Either data output or transformation output files have to be provided!");
      return ILLEGAL_PARAMETERS;
    }
    // check whether number of input files equals number of output files:
    if (!outs.empty() && (ins.size() != outs.size()))
    {
      writeLog_("Error: The number of input and output files has to be equal!");
      return ILLEGAL_PARAMETERS;
    }
    if (!trafos.empty() && (ins.size() != trafos.size()))
    {
      writeLog_("Error: The number of input and transformation output files has to be equal!");
      return ILLEGAL_PARAMETERS;
    }
    // check whether all input files have the same type (this type is used to store the output type too):
    FileTypes::Type in_type = FileHandler::getType(ins[0]);
    for (Size i = 1; i < ins.size(); ++i)
    {
      if (FileHandler::getType(ins[i]) != in_type)
      {
        writeLog_("Error: All input files have to be in the same format!");
        return ILLEGAL_PARAMETERS;
      }
    }

    if (check_ref) // a valid index OR file should be given
    {
      Size reference_index = getIntOption_("reference:index");
      String reference_file = getStringOption_("reference:file");
      if (reference_index > getStringList_("in").size())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'reference:index' must not be higher than the number of input files");
      }
      if (reference_index && !reference_file.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "'reference:index' and 'reference:file' cannot be used together");
      }

      // file should have same type as other input
      if (!reference_file.empty())
      {
        if (FileHandler::getType(reference_file) != in_type)
        {
          writeLog_("Error: Reference file has not the same format as other input files!");
          return ILLEGAL_PARAMETERS;
        }
      }
    }

    //-------------------------------------------------------------
    // set up alignment algorithm
    //-------------------------------------------------------------
    Param alignment_param = getParam_().copy("algorithm:", true);

    writeDebug_("Used alignment parameters", alignment_param, 3);
    alignment->setParameters(alignment_param);
    alignment->setLogType(log_type_);

    return EXECUTION_OK;
  }

  /// deprecated? (not used in PoseClustering... and moved to initialize_() )
  ExitCodes commonMain_(MapAlignmentAlgorithm* alignment)
  {
    ExitCodes ret = initialize_(alignment);
    if (ret != EXECUTION_OK) return ret;

    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);


    StringList ins = getStringList_("in");
    StringList outs = getStringList_("out");
    StringList trafos = getStringList_("trafo_out");
    Param model_params = getParam_().copy("model:", true);
    String model_type = model_params.getValue("type");
    model_params = model_params.copy(model_type + ":", true);
    FileTypes::Type in_type = FileHandler::getType(ins[0]);
    std::vector<TransformationDescription> transformations;

    //-------------------------------------------------------------
    // perform peak alignment
    //-------------------------------------------------------------
    if (in_type == FileTypes::MZML)
    {
      // load input
      std::vector<MSExperiment<> > peak_maps(ins.size());
      MzMLFile f;
      f.setLogType(log_type_);
      for (Size i = 0; i < ins.size(); ++i)
      {
        f.load(ins[i], peak_maps[i]);
      }

      // try to align
      try
      {
        alignment->alignPeakMaps(peak_maps, transformations);
      }
      catch (Exception::NotImplemented&)
      {
        writeLog_("Error: The algorithm '" + alignment->getName() + "' cannot be used for peak data!");
        return INTERNAL_ERROR;
      }
      if (model_type != "none")
      {
        alignment->fitModel(model_type, model_params, transformations);
      }

      // write output
      progresslogger.startProgress(0, outs.size(), "applying RT transformations and writing output files");
      for (Size i = 0; i < outs.size(); ++i)
      {
        progresslogger.setProgress(i);

        MapAlignmentTransformer::transformRetentionTimes(peak_maps[i],
                                                         transformations[i]);
        //annotate output with data processing info
        addDataProcessing_(peak_maps[i],
                           getProcessingInfo_(DataProcessing::ALIGNMENT));

        f.store(outs[i], peak_maps[i]);
      }
      progresslogger.endProgress();
    }
    //-------------------------------------------------------------
    // perform feature alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::FEATUREXML)
    {
      // load input
      std::vector<std::vector<Peak2D> > feat_maps(ins.size());
      FeatureXMLFile f;
      // f.setLogType(log_type_); // TODO
      progresslogger.startProgress(0, ins.size(), "loading input files");
      for (Size i = 0; i < ins.size(); ++i)
      {
        progresslogger.setProgress(i);
        FeatureMap feature_map;
        f.load(ins[i], feature_map);
        feat_maps[i].resize(feature_map.size());

        FeatureMap::const_iterator it = feature_map.begin();
        std::vector<Peak2D>::iterator c_it = feat_maps[i].begin();
        for (; it != feature_map.end(); ++it, ++c_it)
        {
          *c_it = reinterpret_cast<const Peak2D&>(*it);
        }
      }
      progresslogger.endProgress();

      // try to align
      try
      {
        alignment->alignCompactFeatureMaps(feat_maps, transformations);
      }
      catch (Exception::NotImplemented&)
      {
        writeLog_("Error: The algorithm '" + alignment->getName() + 
                  "' cannot be used for feature data!");
        return INTERNAL_ERROR;
      }
      if (model_type != "none")
      {
        alignment->fitModel(model_type, model_params, transformations);
      }

      // write output
      progresslogger.startProgress(0, outs.size(), "applying RT transformations and writing output files");
      for (Size i = 0; i < outs.size(); ++i)
      {
        progresslogger.setProgress(i);

        FeatureMap feature_map;
        f.load(ins[i], feature_map);

        MapAlignmentTransformer::transformRetentionTimes(feature_map,
                                                         transformations[i]);

        //annotate output with data processing info
        addDataProcessing_(feature_map, 
                           getProcessingInfo_(DataProcessing::ALIGNMENT));

        f.store(outs[i], feature_map);
      }
      progresslogger.endProgress();
    }
    //-------------------------------------------------------------
    // perform consensus alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::CONSENSUSXML)
    {
      // load input
      std::vector<ConsensusMap> cons_maps(ins.size());
      ConsensusXMLFile f;
      // f.setLogType(log_type_); // TODO
      progresslogger.startProgress(0, ins.size(), "loading input files");
      for (Size i = 0; i < ins.size(); ++i)
      {
        progresslogger.setProgress(i);
        f.load(ins[i], cons_maps[i]);
      }
      progresslogger.endProgress();

      // try to align
      try
      {
        alignment->alignConsensusMaps(cons_maps, transformations);
      }
      catch (Exception::NotImplemented&)
      {
        writeLog_("Error: The algorithm '" + alignment->getName() +
                  "' cannot be used for consensus feature data!");
        return INTERNAL_ERROR;
      }
      if (model_type != "none")
      {
        alignment->fitModel(model_type, model_params, transformations);
      }

      // write output
      progresslogger.startProgress(0, outs.size(), "applying RT transformations and writing output files");
      for (Size i = 0; i < outs.size(); ++i)
      {
        progresslogger.setProgress(i);

        MapAlignmentTransformer::transformRetentionTimes(cons_maps[i],
                                                         transformations[i]);

        //annotate output with data processing info
        addDataProcessing_(cons_maps[i],
                           getProcessingInfo_(DataProcessing::ALIGNMENT));

        f.store(outs[i], cons_maps[i]);
      }
      progresslogger.endProgress();
    }
    //-------------------------------------------------------------
    // perform peptide alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::IDXML)
    {
      // load input
      std::vector<std::vector<ProteinIdentification> >
        protein_ids_vec(ins.size());
      std::vector<std::vector<PeptideIdentification> >
        peptide_ids_vec(ins.size());

      IdXMLFile f;
      // f.setLogType_(log_type_);

      progresslogger.startProgress(0, ins.size(), "loading input files");
      for (Size i = 0; i < ins.size(); ++i)
      {
        progresslogger.setProgress(i);
        f.load(ins[i], protein_ids_vec[i], peptide_ids_vec[i]);
      }
      progresslogger.endProgress();

      // try to align
      try
      {
        alignment->alignPeptideIdentifications(peptide_ids_vec,
                                               transformations);
      }
      catch (Exception::NotImplemented&)
      {
        writeLog_("Error: The algorithm '" + alignment->getName() +
                  "' cannot be used for peptide data!");
        return INTERNAL_ERROR;
      }
      if (model_type != "none")
      {
        alignment->fitModel(model_type, model_params, transformations);
      }

      // write output
      progresslogger.startProgress(0, outs.size(), "applying RT transformations and writing output files");
      for (Size i = 0; i < outs.size(); ++i)
      {
        progresslogger.setProgress(i);
        MapAlignmentTransformer::transformRetentionTimes(peptide_ids_vec[i],
                                                         transformations[i]);
        f.store(outs[i], protein_ids_vec[i], peptide_ids_vec[i]);
      }
      progresslogger.endProgress();
    }
    else
    {
      // TODO can this really happen? I think it is tested above. Otherwise
      // throw an appropriate exception?
      return ILLEGAL_PARAMETERS;
    }

    if (!trafos.empty())
    {
      for (Size i = 0; i < transformations.size(); ++i)
      {
        TransformationXMLFile().store(trafos[i], transformations[i]);
      }
    }

    return EXECUTION_OK;
  }

};

}

/// @endcond

#endif // OPENMS_APPLICATIONS_MAPALIGNERBASE_H
