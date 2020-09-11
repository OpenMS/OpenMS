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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHProFilterAlgorithm.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHDeconv TOPP_FLASHDeconv
  (Need to be modified)

  @brief  @ref
  @code
  @endcode
  @verbinclude
  @htmlinclude
*/
// We do not want this class to show up in the docu:
// NEED to fill this part later


class TOPPFLASHProFilter :
    public TOPPBase
{
public:
  TOPPFLASHProFilter() :
      TOPPBase("TOPPFLASHProFilter", "tmp",
               false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {

  }

  ExitCodes main_(int, const char **) override
  {
    MSExperiment map;
    MzMLFile mzml;
    String infile = "/Users/kyowonjeong/Google Drive/ProteinFilter/myo_707_ETDReagentTarget_1e+06__deconved.mzml";
    String fasta = "/Users/kyowonjeong/Google Drive/ProteinFilter/uniprot-proteome_yeast_UP000002311_Myo.fasta";

    double elapsed_wall_secs = 0;
    //double elapsed_deconv_cpu_secs = 0, elapsed_deconv_wall_secs = 0;

    mzml.setLogType(log_type_);
    mzml.load(infile, map);
    auto flashpro = FLASHProFilterAlgorithm(fasta);
    int scan = 1;
    for (auto &it : map)
    {
      if (it.getMSLevel() < 2)
      {
        continue;
      }
      if (it.size() <= 0)
      {
        continue;
      }
      //std::cout <<it.size()<<std::endl;
      auto t_start = chrono::high_resolution_clock::now();
      auto scores = flashpro.getScores(it, 0);
      elapsed_wall_secs = chrono::duration<double>(
          chrono::high_resolution_clock::now() - t_start).count();
      std::cout << scan++ << " -- done [took " << elapsed_wall_secs
                << " s (Wall)] --"
                << endl;
      delete[] scores;
      if (scan > 100)
      {
        //break;
      }
    }

    return
        EXECUTION_OK;
  }


  // the actual main function needed to create an executable
};

int main(int argc, const char **argv)
{
  TOPPFLASHProFilter tool;
  return tool.main(argc, argv);
}
