// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Oliver Alka $
// $Authors: Leon Bichmann, Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_NovorAdapter NovorAdapter

    @brief Novoradapter does ...

    This tool can be used for scientific stuff.

    And more scientific applications.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NovorAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NovorAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPNovorAdapter :
  public TOPPBase
{
public:
  TOPPNovorAdapter() :
    TOPPBase("NovorAdapter", "Template for Tool creation", false)
  {

  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {

    // input & output
    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerOutputFile_("out", "<file>", "", "Novor idXML output");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    // parameters
    
    registerStringOption_("enzyme", "<choice>", "Trypsin", "Digestion enzyme - currently only Trypsin is supported ", false);
    setValidStrings_("enzyme", ListUtils::create<String>("Trypsin"));
    
    registerStringOption_("fragmentation", "<Choice>", "CID", "Fragmentation method", false);
    setValidStrings_("fragmentation", ListUtils::create<String>("CID,HCD"))
    registerStringOption("massAnalyzer", "<Choice>" , "Trap", "MassAnalyzer e.g. (Oritrap CID-Trap, CID-FT, HCD-FT; QTof CID-TOF)", false );
    setValidStrings("massAnalyzer", ListUtils::create<String>("Trap,TOF,FT"));

    registerIntOption_("fragmentIonErrorTol", "<num>", 0.5, "Fragmentation error tolerance  (Da)");
    registerIntOption_("precursorErrorTol", "<num>" , 15, "Precursor error tolerance  (ppm or Da)");
    //TODO: Add choice ppm or Da; Has to be used a string with 0.5Da and/or 15ppm (stringoption?)  
    
    //TODO: Add List with variable and fixed Mod
    regsi (Check which modifcations are implemented: http://rapidnovor.com/wiki/Built-in_PTMs)  

    //TODO: Add StringList - see Modifications 
    registerStringOption("forbiddenResidues", "")
    
   // parameter novorFile will not be wrapped here
  

  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

  }

};


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPNovorAdapter tool;
  return tool.main(argc, argv);
}
/// @endcond
