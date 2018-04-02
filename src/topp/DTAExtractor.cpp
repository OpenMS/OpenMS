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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_DTAExtractor DTAExtractor

    @brief Extracts scans of an mzML file to several files in DTA format.
<center>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ DTAExtractor \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
        </tr>
    </table>
</center>
    The retention time, the m/z ratio (for MS level > 1) and the file extension are appended to the output file name.
    You can limit the exported spectra by m/z range, retention time range or MS level.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_DTAExtractor.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_DTAExtractor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDTAExtractor :
  public TOPPBase
{
public:
  TOPPDTAExtractor() :
    TOPPBase("DTAExtractor", "Extracts spectra of an MS run file to several files in DTA format.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerStringOption_("out", "<file>", "", "base name of DTA output files (RT, m/z and extension are appended)");
    registerStringOption_("mz", "[min]:[max]", ":", "m/z range of precursor peaks to extract.\n"
                                                    "This option is ignored for MS level 1", false);
    registerStringOption_("rt", "[min]:[max]", ":", "retention time range of spectra to extract", false);
    registerStringOption_("level", "i[,j]...", "1,2,3", "MS levels to extract", false);
  }

  ExitCodes main_(int, const char**) override
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //ranges
    String mz, rt, tmp;
    double mz_l, mz_u, rt_l, rt_u;
    vector<UInt> levels;
    //initialize ranges
    mz_l = rt_l = -numeric_limits<double>::max();
    mz_u = rt_u = numeric_limits<double>::max();

    rt = getStringOption_("rt");
    mz = getStringOption_("mz");
    String level = getStringOption_("level");

    //convert bounds to numbers
    try
    {
      //rt
      parseRange_(rt, rt_l, rt_u);
      writeDebug_("rt lower/upper bound: " + String(rt_l) + " / " + String(rt_u), 1);

      //mz
      parseRange_(mz, mz_l, mz_u);
      writeDebug_("mz lower/upper bound: " + String(mz_l) + " / " + String(mz_u), 1);

      //levels
      tmp = level;
      if (level.has(',')) //several levels given
      {
        vector<String> tmp2;
        level.split(',', tmp2);
        for (vector<String>::iterator it = tmp2.begin(); it != tmp2.end(); ++it)
        {
          levels.push_back(it->toInt());
        }
      }
      else //one level given
      {
        levels.push_back(level.toInt());
      }

      String tmp3("MS levels: ");
      tmp3 = tmp3 + *(levels.begin());
      for (vector<UInt>::iterator it = ++levels.begin(); it != levels.end(); ++it)
      {
        tmp3 = tmp3 + ", " + *it;
      }
      writeDebug_(tmp3, 1);
    }
    catch (Exception::ConversionError& /*e*/)
    {
      writeLog_(String("Invalid boundary '") + tmp + "' given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    MzMLFile f;
    f.setLogType(log_type_);
    f.getOptions().setRTRange(DRange<1>(rt_l, rt_u));
    f.load(in, exp);

    DTAFile dta;

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    for (PeakMap::iterator it = exp.begin(); it != exp.end(); ++it)
    {
      // check for MS-level
      if (std::find(levels.begin(), levels.end(), it->getMSLevel()) == levels.end())
      {
        continue;
      }

      // store spectra
      if (it->getMSLevel() > 1)
      {
        double mz_value = 0.0;
        if (!it->getPrecursors().empty()) mz_value = it->getPrecursors()[0].getMZ();
        if (mz_value < mz_l || mz_value > mz_u)
        {
          continue;
        }
        dta.store(out + "_RT" + String(it->getRT()) + "_MZ" + String(mz_value) + ".dta", *it);
      }
      else
      {
        dta.store(out + "_RT" + String(it->getRT()) + ".dta", *it);
      }
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPDTAExtractor tool;
  return tool.main(argc, argv);
}

/// @endcond
