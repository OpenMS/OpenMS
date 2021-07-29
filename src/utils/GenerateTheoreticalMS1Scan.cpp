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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <iostream>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_GenerateTheoreticalMS1Scan GenerateTheoreticalMS1Scan

    @brief Generate Multiple MS1 scans from masses and charge range
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class GenerateTheoreticalMS1Scan :
    public TOPPBase
{
public:
  GenerateTheoreticalMS1Scan() :
      TOPPBase("GenerateTheoreticalMS1Scan", "Generate Multiple MS1 scans from masses and charge range.", false, {}, false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
//    registerInputFile_("in", "<file>", "", "Input mass list file");
//    setValidFormats_("in", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output mzML file");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
  }

  MSExperiment generateMS1Scans(double mass, int min_charge, int max_charge, int num_of_scans)
  {
    // output
    MSExperiment exp;

    // experimental settings to avoid error
    std::vector<SourceFile> source_file_vec;
    SourceFile sf;
//    sf.setNativeIDTypeAccession("MS:1000768");
    source_file_vec.push_back(sf);
    exp.setSourceFiles(source_file_vec);

    // isotope pattern
    auto generator = new CoarseIsotopePatternGenerator();
    auto iso = generator->estimateFromPeptideWeight(mass);

    for (int s_index = 0; s_index < num_of_scans; ++s_index)
    {
      MSSpectrum spec;
      spec.setMSLevel(1);
//      spec.setType(OpenMS::SpectrumSettings::CENTROID);
      spec.setRT(rt_starts + rt_intervals * s_index);
      std::string native_id = "spectrum=" + std::to_string(s_index+1);
      spec.setNativeID(native_id);

      for (int cs = min_charge; cs <= max_charge; ++cs)
      {
//        double mono_mz = mass/cs + Constants::PROTON_MASS_U;
        for (auto& i : iso)
        {
          Peak1D tmp_p;
          tmp_p.setMZ(i.getMZ()/cs + + Constants::PROTON_MASS_U);
          tmp_p.setIntensity(i.getIntensity() * 10e7);
          spec.push_back(tmp_p);
        }
      }
      exp.addSpectrum(spec);
    }

    return exp;
  }


  ExitCodes main_(int, const char **) override
  {

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
//    String in = getStringOption_("in");
    String out_path = getStringOption_("out");

    //-------------------------------------------------------------
    // reading input
    //------------------------------------------------------------
    double mass = 18000.0;
    MSExperiment exp;
    exp = generateMS1Scans(mass, 15, 18, 4);

    //-------------------------------------------------------------
    // writing output
    //------------------------------------------------------------
    MzMLFile outfile;
    outfile.store(out_path, exp);

    return EXECUTION_OK;
  }

  double rt_starts = 50.0;
  double rt_intervals = 4;

};

int main(int argc, const char ** argv)
{
  GenerateTheoreticalMS1Scan tool;
  return tool.main(argc, argv);
}

/// @endcond
