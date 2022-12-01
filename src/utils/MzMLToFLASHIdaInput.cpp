// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors:  Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MzMLToFLASHIdaInput
    @brief  convert mzML files to FLASHIda test input format
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPP_MzMLToFLASHIdaInput : public TOPPBase
{
public:
  TOPP_MzMLToFLASHIdaInput() : TOPPBase("MzMLToFLASHIdaInput", "Convert mzML file to FLASHIda test input format file", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output txt file", true);
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in"), out = getStringOption_("out");

    PeakMap experiment;
    MzMLFile().load(in, experiment);
    fstream out_stream;
    out_stream.open(out, fstream::out);


    out_stream<<"RT(min)\tmz\n";
    for(auto& it : experiment)
    {
      if (it.getMSLevel() != 2)
      {
        continue;
      }
      out_stream<< it.getRT()/60<<"\t"<<it.getPrecursors()[0].getMZ()<<"\n";


    }
    /*
    for(auto& it : experiment)
    {
      if(it.empty() || it.getMSLevel() != 1)
      {
        continue;
      }
      out_stream << "Spec\t"<<it.getRT()<<"\n";
      for(auto& p:it)
      {
        if(p.getIntensity() <= 0)
        {
          continue;
        }
        out_stream <<p.getMZ()<<"\t"<<p.getIntensity()<<"\n";
      }
    }
*/
    out_stream.close();

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPP_MzMLToFLASHIdaInput tool;
  return tool.main(argc, argv);
}

/// @endcond
