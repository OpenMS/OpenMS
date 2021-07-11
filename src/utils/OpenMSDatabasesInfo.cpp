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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_OpenMSDatabasesInfo OpenMSDatabasesInfo

    @brief Information about OpenMS' internal databases

    This util prints the content of OpenMS' enzyme and modification databases to a TSV file.
    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_OpenMSDatabasesInfo.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_OpenMSDatabasesInfo.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class OpenMSDatabasesInfo final :
    public TOPPBase
{
public:
  OpenMSDatabasesInfo() :
      TOPPBase("OpenMSDatabasesInfo", "Prints the content of OpenMS' enzyme and modification databases to TSV", false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    // Output CSV file
    registerOutputFile_("enzymes_out", "<out>", "", "Currently supported enzymes as TSV", true, false);
    setValidFormats_("enzymes_out", ListUtils::create<String>("tsv"));
    registerOutputFile_("mods_out", "<out>", "", "Currently supported modifications as TSV", true, false);
    setValidFormats_("mods_out", ListUtils::create<String>("tsv"));
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**) final override
  {
    auto* enz_db = ProteaseDB::getInstance();
    enz_db->writeTSV(getStringOption_("enzymes_out"));

    auto* mod_db = ModificationsDB::getInstance();
    mod_db->writeTSV(getStringOption_("mods_out"));
    
    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char **argv)
{
  OpenMSDatabasesInfo tool;
  return tool.main(argc, argv);
}


/// @endcond
