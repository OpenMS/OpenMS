// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MapNormalizer MapNormalizer

  @brief Normalizes peak intensities to the percentage of the maximum intensity in the HPLC-MS map.
  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapNormalizer \f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
  </tr>
  </table>
  </center>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MapNormalizer.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_MapNormalizer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapNormalizer :
  public TOPPBase
{
public:
  TOPPMapNormalizer() :
    TOPPBase("MapNormalizer", "Normalizes peak intensities in an MS run.")
  {

  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
  }

  ExitCodes main_(int, const char **) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    MzMLFile f;
    f.load(in, exp);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    //determine maximum peak
    exp.updateRanges();
    double max = exp.getMaxInt() / 100.0;

    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      if (it->getMSLevel() < 2)
      {
        for (PeakMap::SpectrumType::Iterator it2 = it->begin(); it2 != it->end(); ++it2)
        {
          it2->setIntensity(it2->getIntensity() / max);
        }
      }
    }


    /// @todo add chromatogram support for normalization, e.g. for MRM stuff (Andreas)
    /*
      vector<MSChromatogram > chroms = exp.getChromatograms();
      double sum(0);
for (vector<MSChromatogram >::iterator it = chroms.begin(); it != chroms.end(); ++it)
{
  for (MSChromatogram::Iterator it2 = it->begin(); it2 != it->end(); ++it2)
  {
              sum += it2->getIntensity();
          }
      }

      for (vector<MSChromatogram >::iterator it = chroms.begin(); it != chroms.end(); ++it)
      {
          for (MSChromatogram::Iterator it2 = it->begin(); it2 != it->end(); ++it2)
          {
              it2->setIntensity(it2->getIntensity() / sum * 1000000.0);
          }
      }

      exp.setChromatograms(chroms);
    */

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::NORMALIZATION));

    f.store(out, exp);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPMapNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
