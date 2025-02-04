// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDRTCalibration IDRTCalibration

@brief Can be used to calibrate the RTs of peptide hits linearly to standards.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; IDRTCalibration &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer (or other tools operating @n with identifications (in idXML format))</td>
        </tr>
    </table>
</CENTER>

This tool can be used to linearly align RTs of the idXML-File to a reference. If only calibrant_1_input and
calibrant_2_input are given, the first calibrant will result at RT 0.1 and calibrant_2_input will be at 0.9.
If one wants to align the RTs of this idXML file to the IDs of a reference file one can also give the RTs
of the same calibrant in the reference file (calibrant_1_reference, calibrant_2_reference). If these calibrants
are given, the linear transformation (shift and scale) will be calculated such that calibrant_1_input will
be at the same RT as calibrant_1_reference and calibrant_2_input will
be at the same RT as calibrant_2_reference. This only applies if calibrant_1* has a smaller RT than calibrant_2*.
Otherwise the values are swapped.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDRTCalibration.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDRTCalibration.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDRTCalibration :
  public TOPPBase
{
public:
  TOPPIDRTCalibration() :
    TOPPBase("IDRTCalibration", "Calibrate RTs of peptide hits linearly to standards.")
  {

  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerDoubleOption_("calibrant_1_reference", "<RT>", 0.1, "The RT of the first calibrant in the reference file.", false);
    registerDoubleOption_("calibrant_2_reference", "<RT>", 0.9, "The RT of the second calibrant in the reference file.", false);
    registerDoubleOption_("calibrant_1_input", "<RT>", -1.0, "The RT of the first calibrant in the input file. Please note that this value needs to be set. The default value -1.0 is not allowed.", false);
    registerDoubleOption_("calibrant_2_input", "<RT>", -1.0, "The RT of the second calibrant in the input file. Please note that this value needs to be set. The default value -1.0 is not allowed.", false);
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in_file = getStringOption_("in");
    String out_file = getStringOption_("out");

    double rt_calibrant_1_input = getDoubleOption_("calibrant_1_input");
    double rt_calibrant_2_input =  getDoubleOption_("calibrant_2_input");
    double rt_calibrant_1_reference =  getDoubleOption_("calibrant_1_reference");
    double rt_calibrant_2_reference =  getDoubleOption_("calibrant_2_reference");

    if (rt_calibrant_1_input == rt_calibrant_2_input)
    {
      OPENMS_LOG_ERROR << "rt_calibrant_1_input and rt_calibrant_2_input must not have the same value";
      return ILLEGAL_PARAMETERS;
    }
    if (rt_calibrant_1_reference == rt_calibrant_2_reference)
    {
      OPENMS_LOG_ERROR << "rt_calibrant_1_reference and rt_calibrant_2_reference must not have the same value";
      return ILLEGAL_PARAMETERS;
    }
    if (rt_calibrant_1_reference == -1 || rt_calibrant_2_reference == -1)
    {
      OPENMS_LOG_ERROR << "rt_calibrant_1_reference and rt_calibrant_2_reference must be set";
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // testing whether input and output files are accessible
    //-------------------------------------------------------------

    if (rt_calibrant_1_input > rt_calibrant_2_input)
    {
      double temp = rt_calibrant_1_input;
      rt_calibrant_1_input = rt_calibrant_2_input;
      rt_calibrant_2_input = temp;
    }
    if (rt_calibrant_1_reference > rt_calibrant_2_reference)
    {
      double temp = rt_calibrant_1_reference;
      rt_calibrant_1_reference = rt_calibrant_2_reference;
      rt_calibrant_2_reference = temp;
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    FileHandler file;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    file.loadIdentifications(in_file, protein_identifications, identifications, {FileTypes::IDXML});

    for (Size i = 0; i < identifications.size(); ++i)
    {
      if (identifications[i].hasRT())
      {
        double temp_rt = identifications[i].getRT();
        temp_rt = (temp_rt - rt_calibrant_1_input) / (rt_calibrant_2_input - rt_calibrant_1_input)
                  * (rt_calibrant_2_reference - rt_calibrant_1_reference) + rt_calibrant_1_reference;
        identifications[i].setRT(temp_rt);
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    file.storeIdentifications(out_file,
               protein_identifications,
               identifications, {FileTypes::IDXML});

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDRTCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
