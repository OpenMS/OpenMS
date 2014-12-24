// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl, Chris Bielow $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Factory.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MapAlignmentEvaluation MapAlignmentEvaluation

    @brief Evaluates alignment results against a ground truth.
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapAlignmentEvaluationAlgorithm \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (text output)</td>
        </tr>
  </table>
</CENTER>

    This tool implements the evaluation measures published in\n
    "Critical assessment of alignment procedures for LC-MS proteomics and metabolomics measurements",\n
    Eva Lange, Ralf Tautenhahn, Steffen Neumann, Clemens Groepl. BMC Bioinformatics 2008, 9:375.\n
    doi:10.1186/1471-2105-9-375.\n

    Input is a ground truth file as described on the CAAP web page\n
    Output is a recall- or a precision-value.\n

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MapAlignmentEvaluation.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MapAlignmentEvaluation.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignmentEvaluation :
  public TOPPBase
{

public:
  TOPPMapAlignmentEvaluation() :
    TOPPBase("MapAlignmentEvaluation", "Evaluates alignment results against a ground truth.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file: tool", true);
    setValidFormats_("in", ListUtils::create<String>("consensusXML"));
    registerInputFile_("gt", "<file>", "", "input file: ground truth", true);
    setValidFormats_("gt", ListUtils::create<String>("consensusXML"));
    registerStringOption_("type", "<name>", "", "Caap Evaluation type", true);
    StringList types = Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts();
    types.push_back("F1");
    setValidStrings_("type", types);
    registerDoubleOption_("rt_dev", "<double>", 0.1, "Maximum allowed deviation of the retention time", false);
    registerDoubleOption_("mz_dev", "<double>", 0.1, "Maximum allowed deviation of m/z", false);
    registerDoubleOption_("int_dev", "<double>", 100, "Maximum allowed deviation of Intensity", false);
    registerFlag_("use_charge", "Use charge criterion when assesing if two features are identical.", false);
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in   = getStringOption_("in");
    String gt   = getStringOption_("gt");
    String type = getStringOption_("type");

    double rt_dev = getDoubleOption_("rt_dev");
    double mz_dev = getDoubleOption_("mz_dev");
    double int_dev = getDoubleOption_("int_dev");

    bool use_charge = getFlag_("use_charge");

    //-------------------------------------------------------------
    // check for valid input
    //-------------------------------------------------------------
    //check if both input files have the correct type
    if (FileHandler::getType(in) != FileTypes::CONSENSUSXML)
    {
      writeLog_("Error: The input file must be of type consensusXML!");
      return ILLEGAL_PARAMETERS;
    }

    if (FileHandler::getType(gt) != FileTypes::CONSENSUSXML)
    {
      writeLog_("Error: The groundtruth file must be of type consensusXML!");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // read input files
    //-------------------------------------------------------------

    // reader
    ConsensusXMLFile consensus_xml_file_in;
    consensus_xml_file_in.setLogType(log_type_);

    // tool -> consensus_map_in
    ConsensusMap consensus_map_in;
    consensus_xml_file_in.load(in, consensus_map_in);

    // gt -> consensus_map_gt
    ConsensusMap consensus_map_gt;
    consensus_xml_file_in.load(gt, consensus_map_gt);


    //-------------------------------------------------------------
    // set up algorithm
    //-------------------------------------------------------------
    if (type == "F1")
    {
      MapAlignmentEvaluationAlgorithm * algorithm_p = Factory<MapAlignmentEvaluationAlgorithm>::create("precision");
      MapAlignmentEvaluationAlgorithm * algorithm_r = Factory<MapAlignmentEvaluationAlgorithm>::create("recall");

      double precision = 0;
      double recall = 0;

      //evaluate
      algorithm_p->evaluate(consensus_map_in, consensus_map_gt, rt_dev, mz_dev, int_dev, use_charge, precision);
      algorithm_r->evaluate(consensus_map_in, consensus_map_gt, rt_dev, mz_dev, int_dev, use_charge, recall);

      //write output
      cout << "precision" << ": " << precision << "\n";
      cout << "   recall" << ": " << recall << "\n";
      cout << "-->    F1" << ": " << (2 * precision * recall) / (precision + recall) << " (2*precision*recall)/(precision+recall)\n";

    }
    else
    {
      MapAlignmentEvaluationAlgorithm * algorithm = Factory<MapAlignmentEvaluationAlgorithm>::create(type);

      double result = 0;

      //evaluate
      algorithm->evaluate(consensus_map_in, consensus_map_gt, rt_dev, mz_dev, int_dev, use_charge, result);

      //write output
      cout << type << ": " << result << "\n";
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPMapAlignmentEvaluation tool;
  return tool.main(argc, argv);
}

/// @endcond
