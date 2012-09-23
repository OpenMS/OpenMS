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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_FuzzyDiff FuzzyDiff

    @brief Compares two files, tolerating numeric differences.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_FuzzyDiff.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_FuzzyDiff.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFuzzyDiff :
  public TOPPBase
{
public:
  TOPPFuzzyDiff() :
    TOPPBase("FuzzyDiff", "Compares two files, tolerating numeric differences.", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    addEmptyLine_();
    addText_("Input files:");
    registerInputFile_("in1", "<file>", "", "first input file", true, false);
    registerInputFile_("in2", "<file>", "", "second input file", true, false);
    addEmptyLine_();
    addText_("Allowed numeric differences:");
    registerDoubleOption_("ratio", "<double>", 1, "acceptable relative error", false, false);
    setMinFloat_("ratio", 1);
    registerDoubleOption_("absdiff", "<double>", 0, "acceptable absolute difference", false, false);
    setMinFloat_("absdiff", 0);
    addText_("Only one of the criteria has to be satisfied.  Use \"absdiff\" to deal with cases like \"zero vs. epsilon\".");
    addEmptyLine_();
    registerStringList_("whitelist", "<string list>", StringList::create("<?xml-stylesheet"), "Lines containing one of these strings are skipped", false, true);
    addEmptyLine_();
    addText_("Output style:");
    registerIntOption_("verbose", "<int>", 2, "set verbose level:\n"
                                              "0 = very quiet mode (absolutely no output)\n"
                                              "1 = quiet mode (no output unless differences detected)\n"
                                              "2 = default (include summary at end)\n"
                                              "3 = continue after errors\n",
                       false, false
                       );
    setMinInt_("verbose", 0);
    setMaxInt_("verbose", 3);
    registerIntOption_("tab_width", "<int>", 8, "tabulator width, used for calculation of column numbers", false, false);
    setMinInt_("tab_width", 1);
    registerIntOption_("first_column", "<int>", 1, "number of first column, used for calculation of column numbers", false, false);
    setMinInt_("first_column", 0);
    addText_("In the diff output, \"position\" refers to the characters in the string, whereas \"column\" is meant for the text editor.");
  }

  ExitCodes main_(int, const char **)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in1 = getStringOption_("in1");
    String in2 = getStringOption_("in2");
    double acceptable_ratio = getDoubleOption_("ratio");
    double acceptable_absdiff = getDoubleOption_("absdiff");
    StringList whitelist = getStringList_("whitelist");
    int verbose_level = getIntOption_("verbose");
    int tab_width = getIntOption_("tab_width");
    int first_column = getIntOption_("first_column");

    // This is for debugging the parsing of whitelist_ from cmdline or ini file.  Converting StringList back to String is intentional.
    writeDebug_(String("whitelist: ") + String(whitelist) + " (size: " + whitelist.size() + ")", 1);

    OpenMS::FuzzyStringComparator fsc;

    fsc.setAcceptableRelative(acceptable_ratio);
    fsc.setAcceptableAbsolute(acceptable_absdiff);
    fsc.setWhitelist(whitelist);
    fsc.setVerboseLevel(verbose_level);
    fsc.setTabWidth(tab_width);
    fsc.setFirstColumn(first_column);

    if (fsc.compareFiles(in1, in2))
    {
      return EXECUTION_OK;
    }
    else
    {
      // TODO think about better exit codes.
      return PARSE_ERROR;
    }
  }

};

int main(int argc, const char ** argv)
{
  TOPPFuzzyDiff tool;
  return tool.main(argc, argv);
}

/// @endcond
