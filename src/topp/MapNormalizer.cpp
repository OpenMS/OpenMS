// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
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
@page TOPP_MapNormalizer MapNormalizer

@brief Normalizes peak intensities to the percentage of the maximum intensity in the HPLC-MS map.
<center>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=2> &rarr; MapNormalizer &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
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
    FileHandler f;
    f.loadExperiment(in, exp, {FileTypes::MZML});

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    //determine maximum peak
    exp.updateRanges();
    double max = exp.getMaxIntensity() / 100.0;

    for (MSSpectrum& ms : exp)
    {
      if (ms.getMSLevel() < 2)
      {
        for (Peak1D& pk : ms)
        {
          pk.setIntensity(pk.getIntensity() / max);
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

    f.storeExperiment(out, exp, {FileTypes::MZML});

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPMapNormalizer tool;
  return tool.main(argc, argv);
}

/// @endcond
