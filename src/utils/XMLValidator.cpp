// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_XMLValidator XMLValidator

    @brief Validates XML files against an XSD schema.

    When a schema file is given, the input file is simply validated against the schema.

    When no schema file is given, the tool tries to determine the file type and
    validates the file against the latest schema version.

    @note XML schema files for the %OpenMS XML formats and several other XML
    formats can be found in the folder
          OpenMS/share/OpenMS/SCHEMAS/

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_XMLValidator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_XMLValidator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPXMLValidator :
  public TOPPBase
{
public:
  TOPPXMLValidator() :
    TOPPBase("XMLValidator", "Validates XML files against an XSD schema.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "file to validate");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzData,featureXML,mzid,idXML,consensusXML,mzXML,ini,pepXML,TraML,xml"));
    registerInputFile_("schema", "<file>", "", "schema to validate against.\nIf no schema is given, the file is validated against the latest schema of the file type.", false);
    setValidFormats_("schema", ListUtils::create<String>("xsd"));
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String schema = getStringOption_("schema");
    bool valid = false;

    if (!schema.empty()) //schema explicitly given
    {
      valid = XMLValidator().isValid(in, schema);
    }
    else //no schema given
    {
      //determine input type
      FileTypes::Type in_type = FileHandler::getType(in);
      if (in_type == FileTypes::UNKNOWN)
      {
        writeLog_("Error: Could not determine input file type and no xsd schema was provided!");
        return PARSE_ERROR;
      }

      cout << endl << "Validating " << FileTypes::typeToName(in_type) << " file";
      switch (in_type)
      {
      case FileTypes::MZDATA:
        cout << " against schema version " << MzDataFile().getVersion() << endl;
        valid = MzDataFile().isValid(in, std::cerr);
        break;

      case FileTypes::FEATUREXML:
        cout << " against schema version " << FeatureXMLFile().getVersion() << endl;
        valid = FeatureXMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::IDXML:
        cout << " against schema version " << IdXMLFile().getVersion() << endl;
        valid = IdXMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::CONSENSUSXML:
        cout << " against schema version " << ConsensusXMLFile().getVersion() << endl;
        valid = ConsensusXMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::MZXML:
        cout << " against schema version " << MzXMLFile().getVersion() << endl;
        valid = MzXMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::INI:
        cout << " against schema version " << ParamXMLFile().getVersion() << endl;
        valid = ParamXMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::PEPXML:
        cout << " against schema version " << PepXMLFile().getVersion() << endl;
        valid = PepXMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::MZML:
        cout << " against schema version " << MzMLFile().getVersion() << endl;
        valid = MzMLFile().isValid(in, std::cerr);
        break;

      case FileTypes::TRAML:
        cout << " against schema version " << TraMLFile().getVersion() << endl;
        valid = TraMLFile().isValid(in, std::cerr);
        break;

      default:
        cout << endl << "Aborted: Validation of this file type is not supported!" << endl;
        return PARSE_ERROR;
      }
    }

    //Result
    if (valid)
    {
      cout << "Success: the file is valid!" << endl;
      return EXECUTION_OK;
    }
    else
    {
      cout << "Failed: errors are listed above!" << endl;
      return PARSE_ERROR;
    }
  }

};

int main(int argc, const char** argv)
{
  TOPPXMLValidator tool;
  return tool.main(argc, argv);
}

/// @endcond
