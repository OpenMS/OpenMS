// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Tom Waschischeck $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <cstdio>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DatabaseSuitability DatabaseSuitability

@brief Calculates a suitability for a database which was used a for peptide identification search. Also reports the quality of LC-MS spectra.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class DatabaseSuitability :
  public TOPPBase
{
public:
  DatabaseSuitability() :
    TOPPBase("DatabaseSuitability", "Computes a suitability score for a database which was used for a peptide identification search. Also reports the quality of LC-MS spectra.", false)
  {
  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_id", "<file>", "", "Input idXML file from peptide search (after FDR)");
    setValidFormats_("in_id", { "idXML" });
    registerInputFile_("in_spec", "<file>", "", "Input MzML file");
    setValidFormats_("in_spec", { "mzML" });
    registerInputFile_("in_novo", "<file>", "", "Input idXML file containing de novo peptides");
    setValidFormats_("in_novo", { "idXML" });
  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in_id = getStringOption_("in_id");
    String in_spec = getStringOption_("in_spec");
    String in_novo = getStringOption_("in_novo");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    IdXMLFile x;
    vector<ProteinIdentification> prot_ids;
    vector<PeptideIdentification> pep_ids;
    x.load(in_id, prot_ids, pep_ids);



    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

  }

};


// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  DatabaseSuitability tool;
  return tool.main(argc, argv);
}
