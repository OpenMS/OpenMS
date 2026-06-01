// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SemanticValidator SemanticValidator

@brief SemanticValidator for XML files which can be semantically validated.

This tool is able to semantically validate an XML file against a CV-mapping
file. The CV-mapping file describes the validity of CV-terms for a given
tag inside the XML. The CV-mapping file conforms to the CvMapping XML
schema (found at /share/OpenMS/SCHEMAS/CvMapping.xsd or
http://www.psidev.info/sites/default/files/CvMapping.xsd). 

Example files that can be semantically validated using this tool are mzML,
TraML, mzIdentML, mzData or any XML file.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SemanticValidator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SemanticValidator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSemanticValidator :
  public TOPPBase
{
public:
  TOPPSemanticValidator() :
    TOPPBase("SemanticValidator", "SemanticValidator for semantically validating certain XML files.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (any xml file)");
    setValidFormats_("in", ListUtils::create<String>("analysisXML,mzML,traML,mzid,mzData,xml"));

    registerInputFile_("mapping_file", "<file>", "", "Mapping file which is used to semantically validate the given XML file against this mapping file (see 'share/OpenMS/MAPPING' for templates).");
    setValidFormats_("mapping_file", ListUtils::create<String>("xml"));

    registerInputFileList_("cv", "<files>", ListUtils::create<String>(""), "Controlled Vocabulary files containg the CV terms (if left empty, a set of default files are used)", false);
    setValidFormats_("cv", ListUtils::create<String>("obo"));
  }

  ExitCodes main_(int, const char**) override
  {
    String in_file = getStringOption_("in");
    String mapping_file = getStringOption_("mapping_file");
    StringList cv_list = getStringList_("cv");

    CVMappings mappings;
    CVMappingFile().load(mapping_file, mappings, false);

    // Allow definition of the controlled vocabulary files on the commandlines.
    // If none are defined, the hardcoded obo files are used
    ControlledVocabulary cv;
    if (!cv_list.empty())
    {
      for (Size i = 0; i < cv_list.size(); i++)
      {
        // TODO do we need to provide the name of the namespace here?
        cv.loadFromOBO("", cv_list[i]);
      }
    }
    else
    {
      cv.loadFromOBO("PSI-MOD", File::find("/CHEMISTRY/PSI-MOD.obo"));
      cv.loadFromOBO("PATO", File::find("/CV/quality.obo"));
      cv.loadFromOBO("UO", File::find("/CV/unit.obo"));
      cv.loadFromOBO("brenda", File::find("/CV/brenda.obo"));
      cv.loadFromOBO("GO", File::find("/CV/goslim_goa.obo"));
      cv.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));
      cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    }

    // check cv params
    Internal::SemanticValidator semantic_validator(mappings, cv);
    semantic_validator.setCheckTermValueTypes(true);
    semantic_validator.setCheckUnits(true);

    StringList errors, warnings;

    bool valid = semantic_validator.validate(in_file, errors, warnings);
    for (Size i = 0; i < warnings.size(); ++i)
    {
      cout << "Warning: " << warnings[i] << endl;
    }
    for (Size i = 0; i < errors.size(); ++i)
    {
      cout << "Error: " << errors[i] << endl;
    }

    if (valid && warnings.empty() && errors.empty())
    {
      cout << "Congratulations, the file is valid!" << endl;
      return EXECUTION_OK;
    }
    else
    {
      return PARSE_ERROR;
    }

  }

};

int main(int argc, const char** argv)
{
  TOPPSemanticValidator tool;
  return tool.main(argc, argv);
}

/// @endcond
