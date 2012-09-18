// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MzTabExporter MzTabExporter

   @brief This application converts several %OpenMS XML formats (featureXML, consensusXML, and idXML) to mzTab.

  <CENTER>
    <table>
     <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MzTabExporter \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> ProteinQuantifier </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> external tools (MS Excel, OpenOffice, Notepad)</td>
    </tr>
   </table>
  </CENTER>

  See the mzTab specification for details on the format.

  @experimental This algorithm and underlying format is work in progress and might change.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MzTabExporter.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_MzTabExporter.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
class TOPPMzTabExporter : public TOPPBase
{
public:
  TOPPMzTabExporter() :
    TOPPBase("MzTabExporter", "Exports various XML formats to an mzTab file.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file annotated by ProteinQuantifier");
    setValidFormats_("in", StringList::create("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
  } 


  ExitCodes main_( int, const char** )
  {
    // parameter handling
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    // input file type
    FileTypes::Type in_type = FileHandler::getType(in);
    writeDebug_(String("Input file type: ") + FileHandler::typeToName(in_type), 2);

    if ( in_type == FileTypes::UNKNOWN )
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    String document_id;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    IdXMLFile().load( in, prot_ids, pep_ids, document_id);

    MzTabFile mztab;
    mztab.store(out, prot_ids, pep_ids, in, document_id);

    return EXECUTION_OK;
  }
};
}

int main( int argc, const char** argv )
{
  TOPPMzTabExporter t;
  return t.main(argc, argv);
}

/// @endcond
