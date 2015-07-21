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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
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
  void handleReference_(MapAlignmentAlgorithmIdentification& alignment)
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
    alignment.setReference(reference_index, reference_file);
  }

  ExitCodes checkParameters_(bool check_ref = false)
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
      writeLog_("Error: Either data output or transformation output files have to be provided (parameters 'out'/'trafo_out')");
      return ILLEGAL_PARAMETERS;
    }
    // check whether number of input files equals number of output files:
    if (!outs.empty() && (ins.size() != outs.size()))
    {
      writeLog_("Error: The number of data input and output files has to be equal (parameters 'in'/'out')");
      return ILLEGAL_PARAMETERS;
    }
    if (!trafos.empty() && (ins.size() != trafos.size()))
    {
      writeLog_("Error: The number of data input and transformation output files has to be equal (parameters 'in'/'trafo_out')");
      return ILLEGAL_PARAMETERS;
    }
    // check whether all input files have the same type (this type is used to store the output type too):
    FileTypes::Type in_type = FileHandler::getType(ins[0]);
    for (Size i = 1; i < ins.size(); ++i)
    {
      if (FileHandler::getType(ins[i]) != in_type)
      {
        writeLog_("Error: All input files (parameter 'in') must have the same format!");
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
          writeLog_("Error: Reference file must have the same format as other input files (parameters 'reference:file'/'in')");
          return ILLEGAL_PARAMETERS;
        }
      }
    }

    return EXECUTION_OK;
  }

};

}

/// @endcond

#endif // OPENMS_APPLICATIONS_MAPALIGNERBASE_H
