// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <iomanip>
#include <cctype>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page
*/

/// @cond

class TOPPParseMsAlign :
    public TOPPBase
{
public:
  TOPPParseMsAlign() :
      TOPPBase("ParseMsAlign", ".", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    //registerOutputFileList_("out_spec", "<file for MS1, file for MS2, ...>", {""},
    //                            "output files (tsv) - spectrum level deconvoluted masses per ms level", false);
    //

    registerInputFile_("in", "<file>", {}, "Input file");
    registerIntOption_("del", "delimiter", 0, "0: space 1: tap", false);
    registerOutputFile_("out", "<file>", "", "");
    //setValidFormats_("i1n", ListUtils::create<String>("tsv"));
    //setValidFormats_("i1n", ListUtils::create<String>("tsv"));
  }

  ExitCodes main_(int, const char **) override
  {
    auto infile = getStringOption_("in");
    auto outfile = getStringOption_("out");
    char delim = getIntOption_("del") == 0 ? ' ' : '\t';
    fstream outstream;
    outstream.open(outfile, fstream::out); //
    std::ifstream instream(infile);
    outstream << "v=[";
    String line;
    int pcharge = 0;
    double pmass = 0;
    while (std::getline(instream, line))
    {
      if (line.hasPrefix("PRECURSOR_CHARGE"))
      {
        pcharge = stoi(line.substr(17));
      }

      if (line.hasPrefix("PRECURSOR_MASS"))
      {
        pmass = stod(line.substr(15));
      }
      if (!isdigit(line[0]))
      {
        continue;
      }
      vector<String> tokens;
      std::stringstream tmp_stream(line);
      String str;

      while (getline(tmp_stream, str, delim))
      {
        tokens.push_back(str);
      }
      double mass = stod(tokens[0]);
      int charge = stoi(tokens[2]);
      if (mass > pmass || charge > pcharge)
      {
        continue;
      }
      outstream << mass << "," << charge << "\n";
    }
    outstream << "];\n";
    outstream.close();
    instream.close();
    return EXECUTION_OK;
  }
};

int main(int argc, const char **argv)
{
  TOPPParseMsAlign tool;
  return tool.main(argc, argv);
}

/// @endcond
