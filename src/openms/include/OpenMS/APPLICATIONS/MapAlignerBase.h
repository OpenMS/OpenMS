// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
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
  // @brief stores model defaults for map aligner algorithms
  struct MapAlignerBase
  {
    // "public" so it can be used in DefaultParamHandlerDocumenter to get docu
    static Param getModelDefaults(const String& default_model)
    {
      Param params;
      params.setValue("type", default_model, "Type of model");
      // TODO: avoid referring to each TransformationModel subclass explicitly
      std::vector<std::string> model_types = {"linear","b_spline","lowess","interpolated"};
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

      TransformationModelLowess::getDefaultParameters(model_params);
      params.insert("lowess:", model_params);
      params.setSectionDescription("lowess", "Parameters for 'lowess' model");

      TransformationModelInterpolated::getDefaultParameters(model_params);
      params.insert("interpolated:", model_params);
      params.setSectionDescription("interpolated",
                                  "Parameters for 'interpolated' model");
      return params;
    }
  };


class TOPPMapAlignerBase :
  public TOPPBase, public MapAlignerBase
{

public:
  TOPPMapAlignerBase(String name, String description, bool official = true) :
    TOPPBase(name, description, official), ref_params_(REF_NONE)
  {
  }

protected:

  // Kind of reference parameters that the tool offers:
  // - REF_NONE: no reference
  // - REF_RESTRICTED: reference file must have same type as input files
  // - REF_FLEXIBLE: reference file can have any supported file type
  enum ReferenceParameterKind { REF_NONE, REF_RESTRICTED, REF_FLEXIBLE }
    ref_params_;

  void registerOptionsAndFlagsMapAligners_(const String& file_formats,
                                           enum ReferenceParameterKind ref_params)
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files to align (all must have the same file type)", true);
    setValidFormats_("in", ListUtils::create<String>(file_formats));
    registerOutputFileList_("out", "<files>", StringList(), "Output files (same file type as 'in'). This option or 'trafo_out' has to be provided; they can be used together.", false);
    setValidFormats_("out", ListUtils::create<String>(file_formats));
    registerOutputFileList_("trafo_out", "<files>", StringList(), "Transformation output files. This option or 'out' has to be provided; they can be used together.", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));

    if (ref_params != REF_NONE)
    {
      registerTOPPSubsection_("reference", "Options to define a reference file (use either 'file' or 'index', not both)");
      String description = "File to use as reference";
      if (ref_params == REF_RESTRICTED)
      {
        description += " (same file format as input files required)";
      }
      registerInputFile_("reference:file", "<file>", "",  description, false);
      setValidFormats_("reference:file", ListUtils::create<String>(file_formats));
      registerIntOption_("reference:index", "<number>", 0, "Use one of the input files as reference ('1' for the first file, etc.).\nIf '0', no explicit reference is set - the algorithm will select a reference.", false);
      setMinInt_("reference:index", 0);
    }
    ref_params_ = ref_params;
  }

  ExitCodes checkParameters_()
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
      writeLogError_("Error: Data output or transformation output files have to be provided (parameters 'out'/'trafo_out')");
      return ILLEGAL_PARAMETERS;
    }
    // check whether number of input files equals number of output files:
    if (!outs.empty() && (ins.size() != outs.size()))
    {
      writeLogError_("Error: The number of data input and output files has to be equal (parameters 'in'/'out')");
      return ILLEGAL_PARAMETERS;
    }
    if (!trafos.empty() && (ins.size() != trafos.size()))
    {
      writeLogError_("Error: The number of data input and transformation output files has to be equal (parameters 'in'/'trafo_out')");
      return ILLEGAL_PARAMETERS;
    }
    // check whether all input files have the same type (this type is used to store the output type too):
    FileTypes::Type in_type = FileHandler::getType(ins[0]);
    for (Size i = 1; i < ins.size(); ++i)
    {
      if (FileHandler::getType(ins[i]) != in_type)
      {
        writeLogError_("Error: All input files (parameter 'in') must have the same format!");
        return ILLEGAL_PARAMETERS;
      }
    }

    if (ref_params_ != REF_NONE) // a valid ref. index OR file should be given
    {
      Size reference_index = getIntOption_("reference:index");
      String reference_file = getStringOption_("reference:file");
      if (reference_index > ins.size())
      {
        writeLogError_("Error: Value of parameter 'reference:index' must not be higher than the number of input files");
        return ILLEGAL_PARAMETERS;
      }
      if (reference_index && !reference_file.empty())
      {
        writeLogError_("Error: Parameters 'reference:index' and 'reference:file' cannot be used together");
        return ILLEGAL_PARAMETERS;
      }

      if ((ref_params_ == REF_RESTRICTED) && !reference_file.empty() &&
          (FileHandler::getType(reference_file) != in_type))
      {
        writeLogError_("Error: Reference file must have the same format as other input files (parameters 'reference:file'/'in')");
        return ILLEGAL_PARAMETERS;
      }
    }

    return EXECUTION_OK;
  }

};

}

/// @endcond

