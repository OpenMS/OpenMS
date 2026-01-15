// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <vector>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ExternalCalibration ExternalCalibration

@brief Performs an mass recalibration on an MS experiment using an external calibration function.

<CENTER>
<table>
<tr>
<th ALIGN = "center"> pot. predecessor tools </td>
<td VALIGN="middle" ROWSPAN=3> &rarr; ExternalCalibration &rarr;</td>
<th ALIGN = "center"> pot. successor tools </td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> any tool operating on MS peak data @n (in mzML format)</td>
</tr>
<tr>
<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
</tr>
</table>
</CENTER>

Recalibrates an MS experiment globally using a constant, linear or quadratic calibration on the observed ppm error,
i.e. using offset=-5, slope=0, power=0 assumes the observed data has -5 ppm decalibration, i.e. the observed m/z is too small and should be increased by 5 ppm! Slope and power
are coefficients for the observed m/z, i.e. y = offset + x * slope + x * x * power, where x is the observed m/z and y 
is the resulting mass correction in ppm. Usually slope and offset are very small (< 0.01).
If you only want a 'rough' recalibration, using offset is usually sufficient.

The calibration function is applied to all scans (with the desired level, see below), i.e. time dependent recalibration cannot be modeled.
  
The user can select what MS levels are subjected to calibration.
Spectra with other MS levels remain unchanged.
Calibration must be done once for each mass analyzer.

Either raw or centroided data can be used as input.
<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ExternalCalibration.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ExternalCalibration.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES



class TOPPExternalCalibration :
  public TOPPBase
{
public:
  TOPPExternalCalibration() :
    TOPPBase("ExternalCalibration", "Applies an external mass recalibration.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    // data
    registerInputFile_("in", "<file>", "", "Input peak file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file ");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
        
    addEmptyLine_();

    registerDoubleOption_("offset", "", 0.0, "Mass offset in ppm", false);
    registerDoubleOption_("slope", "", 0.0, "Slope (dependent on m/z)", false);
    registerDoubleOption_("power", "", 0.0, "Power (dependent on m/z)", false);
    
    addEmptyLine_();
    
    registerIntList_("ms_level", "i j ...", ListUtils::create<int>("1,2,3"), "Target MS levels to apply the transformation onto. Scans with other levels remain unchanged.", false);
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out"); 
    
    IntList ms_level = getIntList_("ms_level");

    double offset = getDoubleOption_("offset");
    double slope = getDoubleOption_("slope");
    double power = getDoubleOption_("power");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    // Raw data
    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    MZTrafoModel tm;
    tm.setCoefficients(offset, slope, power);

    InternalCalibration ic;
    ic.setLogType(log_type_);
    ic.applyTransformation(exp, ms_level, tm);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(exp, getProcessingInfo_(DataProcessing::CALIBRATION));

    FileHandler().storeExperiment(out, exp, {FileTypes::MZML}, log_type_);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPExternalCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
