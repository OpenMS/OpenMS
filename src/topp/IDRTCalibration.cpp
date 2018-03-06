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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

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
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ IDRTCalibration \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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
    TOPPBase("IDRTCalibration", "Can be used to calibrate RTs of peptide hits linearly to standards.")
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
      LOG_ERROR << "rt_calibrant_1_input and rt_calibrant_2_input must not have the same value";
      return ILLEGAL_PARAMETERS;
    }
    if (rt_calibrant_1_reference == rt_calibrant_2_reference)
    {
      LOG_ERROR << "rt_calibrant_1_reference and rt_calibrant_2_reference must not have the same value";
      return ILLEGAL_PARAMETERS;
    }
    if (rt_calibrant_1_reference == -1 || rt_calibrant_2_reference == -1)
    {
      LOG_ERROR << "rt_calibrant_1_reference and rt_calibrant_2_reference must be set";
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
    IdXMLFile file;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    String document_id;
    file.load(in_file, protein_identifications, identifications, document_id);

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

    file.store(out_file,
               protein_identifications,
               identifications);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDRTCalibration tool;
  return tool.main(argc, argv);
}

/// @endcond
