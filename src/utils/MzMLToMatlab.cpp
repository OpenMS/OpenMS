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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QFile>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MzMLToMatlab MzMLToMatlab

    @brief Make spectrum variable txt for MATLAB from mzml file

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MzMLToMatlab.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MzMLToMatlab.html
*/

/// @cond TOPPCLASSES

class TOPPMzMLToMatlab :
    public TOPPBase
{
public:
  TOPPMzMLToMatlab() :
      TOPPBase("MzMLToMatlab", "Converts spectra in a mzML file into MATLAB variable txt", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input msalign file");
      registerInputFile_("log", "<file>", "", "Input log file");
      registerOutputFile_("out", "<file>", "", "out msalign file");
    //setValidFormats_("in", ListUtils::create<String>("mzML"));
  }

    std::vector<std::string> tokenise(const std::string &str){
        std::vector<std::string> tokens;
        Size first = 0;
        while(first<str.size()){
            Size second = str.find_first_of(',',first);
            //first has index of start of token
            //second has index of end of token + 1;
            if(second==std::string::npos){
                second = str.size();
            }
            std::string token = str.substr(first, second-first);
            tokens.push_back(token);
            first = second + 1;
        }
        return tokens;
    }

    ExitCodes main_(int, const char **) override
  {
      String in = getStringOption_("in");
      String log = getStringOption_("log");
      String out = getStringOption_("out");

      std::ifstream ins(in);
      std::ifstream logs(log);
      fstream outs;
      outs.open(out, fstream::out);

      std::map<int,std::vector<double>> maps;
      String line;

      while (std::getline(logs, line)) {
        auto tokens = tokenise(line);
        int scan = std::stoi(tokens[0]);
        double mass = std::stod(tokens[1]);
        double charge = std::stod(tokens[2]);
        double intenisty = std::stod(tokens[3]);
        vector<double> val;
        val.push_back(mass);
          val.push_back(charge);
          val.push_back(intenisty);
        maps[scan] = val;
      }

      int scan = 0;
      vector<double> val;
      bool exist = false;
      outs<<"#FLASHIda precursor information used\n";
      while (std::getline(ins, line)) {
          if(line.hasPrefix("SCANS")){
              scan = std::stoi(line.substr(6));
              if(maps.find(scan) != maps.end()){
                val = maps[scan];
                exist = true;
              }else{
                  exist = false;
              }
          }

          if(exist){
            if(line.hasPrefix("PRECURSOR_CHARGE")){
              outs<< "PRECURSOR_CHARGE="<< (int)val[1]<<"\n";
              continue;
            }
            if(line.hasPrefix("PRECURSOR_MASS")){
                outs<< "PRECURSOR_MASS="<< val[0]<<"\n";
                continue;
            }

            if(line.hasPrefix("PRECURSOR_INTENSITY")){
                outs<< "PRECURSOR_INTENSITY="<< val[2]<<"\n";
                continue;
            }
          }
          outs<<line<<"\n";
      }

      logs.close();
      outs.close();
      ins.close();

    /*String out = FileHandler::stripExtension(in);

    fstream matlabOut;
    matlabOut.open(out + ".m", fstream::out);

    MSExperiment map;
    MzMLFile mzml;
    mzml.load(in, map);
    int cntr = 1;
    for (auto &it : map)
    {
      matlabOut << "s" << it.getMSLevel() << "n" << cntr++ << "=[";
      for (auto &p : it)
      {
        matlabOut << std::to_string(p.getMZ()) << " " << std::to_string(p.getIntensity()) << ";";
      }
      matlabOut << "];\n";
    }
    matlabOut.close();*/
    return EXECUTION_OK;
  }
};

int main(int argc, const char **argv)
{
  TOPPMzMLToMatlab tool;
  return tool.main(argc, argv);
}

/// @endcond
