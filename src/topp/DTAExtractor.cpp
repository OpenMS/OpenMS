// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

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
        <th ALIGN = "center"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=2> &rarr; DTAExtractor &rarr;</td>
        <th ALIGN = "center"> pot. successor tools </td>
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
      writeLogError_(String("Invalid boundary '") + tmp + "' given. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    FileHandler f;
    f.getOptions().setRTRange(DRange<1>(rt_l, rt_u));
    f.loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    FileHandler dta;

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    for (MSSpectrum& spec : exp)
    {
      // check for MS-level
      if (std::find(levels.begin(), levels.end(), spec.getMSLevel()) == levels.end())
      {
        continue;
      }

      // store spectra
      if (spec.getMSLevel() > 1)
      {
        double mz_value = 0.0;
        if (!spec.getPrecursors().empty())
        {
          mz_value = spec.getPrecursors()[0].getMZ();
        }
        if (mz_value < mz_l || mz_value > mz_u)
        {
          continue;
        }
        MSExperiment exp;
        exp.addSpectrum(spec);
        dta.storeExperiment(out + "_RT" + String(spec.getRT()) + "_MZ" + String(mz_value) + ".dta", exp);
      }
      else
      {
        MSExperiment exp;
        exp.addSpectrum(spec);
        dta.storeExperiment(out + "_RT" + String(spec.getRT()) + ".dta", exp);
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
