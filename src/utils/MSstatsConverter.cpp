// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/regex.hpp>
#include <OpenMS/FORMAT/MSstatsFile.h>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MSstatsConverter

    @brief Converter to input for MSstats

    This util consumes an ID-mapped consensusXML file and OpenMS experimental design in TSV format to create a CSV file which can subsequently be used as input for the R package MSstats [1].

    [1] M. Choi et al. MSstats: an R package for statistical analysis for quantitative mass spectrometry-based proteomic experiments. Bioinformatics (2014), 30 (17): 2524-2526

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MSstats.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MSstats.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSstatsConverter final :
  public TOPPBase
{
public:

  TOPPMSstatsConverter() :
          TOPPBase("MSstatsConverter", "Converter to input for MSstats", false)
  {

  }

protected:

    // this function will be used to register the tool parameters
    // it gets automatically called on tool execution
    void registerOptionsAndFlags_() final override
    {
      // Input consensusXML
      registerInputFile_(param_in, "<in>", "", "Input consensusXML with peptide intensities",
                         true, false);
      setValidFormats_(param_in, ListUtils::create<String>("consensusXML"), true);

      registerInputFile_(param_in_design, "<in_design>", "", "Experimental Design file", true,
                         false);
      setValidFormats_(param_in_design, ListUtils::create<String>("tsv"), true);

      registerStringOption_(param_method, "<method>",
                            "LFQ",
                            "Method used in the experiment(label free [LFQ], isobaric labeling [ISO]))", false,
                            false);
      setValidStrings_(param_method,
                       ListUtils::create<String>("LFQ,ISO"));

      registerStringOption_(param_msstats_bioreplicate, "<msstats_bioreplicate>",
                            "MSstats_BioReplicate",
                            "Which column in the condition table should be used for MSstats 'BioReplicate'", false,
                            false);
      registerStringOption_(param_msstats_condition, "<msstats_condition>", "MSstats_Condition",
                            "Which column in the condition table should be used for MSstats 'Condition'", false, false);

      registerStringOption_(param_msstats_mixture, "msstats_mixture", "MSstats_Mixture",
                            "Which column in the condition table should be used for MSstats 'Mixture'", false, false);

      // advanced option to overwrite MS file annotations in consensusXML
      registerInputFileList_(param_reannotate_filenames, "<file(s)>", StringList(),
                             "Overwrite MS file names in consensusXML", false, true);

      // Isotope label type
      registerFlag_(param_labeled_reference_peptides,
                    "If set, IsotopeLabelType is 'H', else 'L'");

      // Specifies how peptide ions eluding at different retention times should be resolved
      registerStringOption_(param_retention_time_summarization_method,
                            "<retention_time_summarization_method>", "max",
                            "How undistinguishable peptides at different retention times should be treated", false,
                            true);
      setValidStrings_(param_retention_time_summarization_method,
                       ListUtils::create<String>("manual,max,min,mean,sum"));

      // Output CSV file
      registerOutputFile_(param_out, "<out>", "", "Input CSV file for MSstats.", true, false);
      setValidFormats_(param_out, ListUtils::create<String>("csv"));
    }

    // the main_ function is called after all parameters are read
    ExitCodes main_(int, const char **) final override
    {
      try
      {
        // Input file, must be consensusXML
        const String arg_in(getStringOption_(param_in));
        const FileTypes::Type in_type(FileHandler::getType(arg_in));

        fatalErrorIf_(
                in_type != FileTypes::CONSENSUSXML,
                "Input type is not consensusXML!",
                ILLEGAL_PARAMETERS);
        // Tool arguments
        const String arg_method = getStringOption_(param_method);
        const String arg_out = getStringOption_(param_out);

        // Experimental Design file
        const String arg_in_design = getStringOption_(param_in_design);
        const ExperimentalDesign design = ExperimentalDesignFile::load(arg_in_design, false);
        ExperimentalDesign::SampleSection sampleSection = design.getSampleSection();

        ConsensusMap consensus_map;
        ConsensusXMLFile().load(arg_in, consensus_map);

        StringList reannotate_filenames = getStringList_(param_reannotate_filenames);
        bool is_isotope_label_type = getFlag_(param_labeled_reference_peptides);
        String bioreplicate = getStringOption_(param_msstats_bioreplicate);
        String condition = getStringOption_(param_msstats_condition);
        String mixture = getStringOption_(param_msstats_mixture);
        String retention_time_summarization_method = getStringOption_(param_retention_time_summarization_method);

        MSstatsFile msStatsFile;


        if (arg_method == "LFQ")
        {
            msStatsFile.storeLFQ(arg_out, consensus_map, design,
                                 reannotate_filenames, is_isotope_label_type,
                                 bioreplicate, condition, retention_time_summarization_method);
        }
        else if (arg_method == "ISO")
        {
            msStatsFile.storeISO(arg_out, consensus_map, design,
                                 reannotate_filenames, bioreplicate, condition, 
                                 mixture, retention_time_summarization_method);
        }
        return EXECUTION_OK;
      }
      catch (const ExitCodes &exit_code)
      {
        return exit_code;
      }

    }

    static const String param_in;
    static const String param_in_design;
    static const String param_method;
    static const String param_msstats_bioreplicate;
    static const String param_msstats_condition;
    static const String param_msstats_mixture;
    static const String param_out;
    static const String param_labeled_reference_peptides;
    static const String param_retention_time_summarization_method;
    static const String param_reannotate_filenames;

private:

    static void fatalErrorIf_(const bool error_condition, const String &message, const int exit_code)
    {
      if (error_condition)
      {
        LOG_FATAL_ERROR << "FATAL: " << message << std::endl;
        throw exit_code;
      }
    }
};

const String TOPPMSstatsConverter::param_in = "in";
const String TOPPMSstatsConverter::param_in_design = "in_design";
const String TOPPMSstatsConverter::param_method = "method";
const String TOPPMSstatsConverter::param_msstats_bioreplicate = "msstats_bioreplicate";
const String TOPPMSstatsConverter::param_msstats_condition = "msstats_condition";
const String TOPPMSstatsConverter::param_msstats_mixture = "msstats_mixture";
const String TOPPMSstatsConverter::param_out = "out";
const String TOPPMSstatsConverter::param_labeled_reference_peptides = "labeled_reference_peptides";
const String TOPPMSstatsConverter::param_retention_time_summarization_method = "retention_time_summarization_method";
const String TOPPMSstatsConverter::param_reannotate_filenames = "reannotate_filenames";

// the actual main function needed to create an executable
int main(int argc, const char **argv)
{
  TOPPMSstatsConverter tool;
  return tool.main(argc, argv);
}


/// @endcond
