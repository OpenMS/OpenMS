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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerIterative.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page UTILS_PeakPickerIterative PeakPickerIterative

        @brief A tool for peak detection in profile data. Executes the peak picking with @ref OpenMS::PeakPickerIterative "high_res" algorithm.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ PeakPickerIterative \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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
  @verbinclude UTILS_PeakPickerIterative.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_PeakPickerIterative.html

*/

typedef PeakMap MapType;

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeakPickerIterative
  : public TOPPBase
{
 public:

  TOPPPeakPickerIterative()
    : TOPPBase("PeakPickerIterative","Finds mass spectrometric peaks in profile mass spectra.", false)
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

    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in,exp);
    PeakPickerIterative pp;
    pp.setParameters(picker_param);
    pp.setLogType(log_type_);
    pp.pickExperiment(exp, out_exp);

    addDataProcessing_(out_exp, getProcessingInfo_(DataProcessing::PEAK_PICKING));
    f.store(out,out_exp); 

    return EXECUTION_OK;
  }

};

int main( int argc, const char** argv )
{
  TOPPPeakPickerIterative tool;
  return tool.main(argc,argv);
}

