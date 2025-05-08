// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIterative.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeakPickerIterative PeakPickerIterative
@brief A tool for peak detection in profile data. Executes the peak picking with @ref OpenMS::PeakPickerIterative "high_res" algorithm.

<CENTER>
	<table>
		<tr>
			<th ALIGN = "center"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=4> &rarr; PeakPickerIterative &rarr;</td>
			<th ALIGN = "center"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_BaselineFilter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> any tool operating on MS peak data @n (in mzML format)</td>
		</tr>
		<tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterGaussian </td>
		</tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_NoiseFilterSGolay </td>
    </tr>
	</table>
</CENTER>

The conversion of the ''raw'' ion count data acquired
by the machine into peak lists for further processing
is usually called peak picking. The choice of the algorithm
should mainly depend on the resolution of the data.
As the name implies, the @ref OpenMS::PeakPickerIterative "high_res"
algorithm is fit for high resolution data.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PeakPickerIterative.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PeakPickerIterative.html

*/

typedef PeakMap MapType;

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPickerIterative
  : public TOPPBase
{
 public:

  TOPPPeakPickerIterative()
    : TOPPBase("PeakPickerIterative","Finds mass spectrometric peaks in profile mass spectra.")
  {
  }

 protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in","<file>","","input file ");
    setValidFormats_("in",ListUtils::create<String>("mzML"));

    registerOutputFile_("out","<file>","","output file");
    setValidFormats_("out",ListUtils::create<String>("mzML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String &) const override
  {
    return PeakPickerIterative().getDefaults();
  }

  ExitCodes main_(int , const char**) override
  {

    String in = getStringOption_("in");
    String out = getStringOption_("out");

    MapType exp;
    MapType out_exp;

    Param picker_param = getParam_().copy("algorithm:", true);

    FileHandler().loadExperiment(in,exp, {FileTypes::MZML}, log_type_);
    PeakPickerIterative pp;
    pp.setParameters(picker_param);
    pp.setLogType(log_type_);
    pp.pickExperiment(exp, out_exp);

    addDataProcessing_(out_exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
    FileHandler().storeExperiment(out,out_exp, {FileTypes::MZML}, log_type_); 

    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{
  TOPPPeakPickerIterative tool;
  return tool.main(argc,argv);
}

///@endcond
