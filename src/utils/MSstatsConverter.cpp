// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MSstatsConverter

    @brief Converter to input for MSstats

    This tool consumes (TODO .... ) to create a file which can subsequently be used as an input for MSstats [1].


    [1] M. Choi et al. “MSstats: an R package for statistical analysis for quantitative mass spectrometry-based proteomic experiments.” Bioinformatics (2014), 30 (17): 2524-2526

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MSstats.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MSstats.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSstatsConverter :
  public TOPPBase
{
public:

    static const String param_input_quantitation;
    static const String param_input_identification;

  TOPPMSstatsConverter() :
    TOPPBase("MSstatsConverter", "Converter to input for MSstats", false)
  {

  }

protected:

  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    // Quantification data
    this->registerInputFile_(TOPPMSstatsConverter::param_input_quantitation, "<input_quantitation>", "", "Quantification data", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_input_quantitation, ListUtils::create<String>("consensusXML"));

    // Identification data
    this->registerInputFile_(TOPPMSstatsConverter::param_input_identification, "<input_identification>", "", "Identification", true, false);
    this->setValidFormats_(TOPPMSstatsConverter::param_input_identification, ListUtils::create<String>("idXML"));
  }


  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
      cout << "Hello World!" << endl;

      return EXECUTION_OK;
  }

};

const String TOPPMSstatsConverter::param_input_quantitation = "input_quantitation";
const String TOPPMSstatsConverter::param_input_identification = "input_identification";


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}
/// @endcond
