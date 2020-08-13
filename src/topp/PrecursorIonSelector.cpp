// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelection.h>
#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_PrecursorIonSelector PrecursorIonSelector

    @brief A tool for precursor ion selection based on identification results.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ PrecursorIonSelector \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
        </tr>
    </table>
</CENTER>

    This tool provides a precursor ion selection based on previous MS/MS
    identifications.

    Different strategies can be chosen:
        <table>
        <tr>
            <td><b>DEX</b></td>
            <td>Dynamic exclusion of features with m/z matching predicted tryptic peptides masses of already identified proteins.</td>
        </tr>
        <tr>
          <td><b>SPS</b></td>
            <td>Selection based on score reflecting the feature's suitability for fragmentation.</td>
        </tr>
        <tr>
            <td><b>Downshift</b></td>
            <td>Similar to DEX, but features are not excluded, only ranked down in the feature list</td>
        </tr>
        <tr>
            <td><b>Upshift</b></td>
            <td>Features with m/z matching predicted tryptic peptide masses of unidentified proteins are ranked up.</td>
        </tr>
        <tr>
            <td><b>IPS</b></td>
            <td>Combination of Down- and Upshift.</td>
        </tr>
        <tr>
            <td><b>ILP_IPS</b></td>
            <td>Iterative precursor ion selection using LP formulations.</td>
        </tr>
    </table>

    This method is described in: Zerck, A.  and Nordhoff, E.  and Resemann, A.  and Mirgorodskaya, E.  and Suckau, D.  and Reinert, K.  and Lehrach, H.  and Gobom, J.:
  An iterative strategy for precursor ion selection for LC-MS/MS based shotgun proteomics, J Prot Res, 2009, 8 (7), 3239-3251.

    Given the feature map of the LC-MS run and the identification results
    the tool determines the next precursors. The precursors are ranked
    depending on the chosen strategy.

    It is also possible run a simulation of selection strategies
    on a complete LC-MS/MS run, e.g. to determine what would have been
    the most efficient strategy.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_PrecursorIonSelector.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_PrecursorIonSelector.html

    For the parameters of the algorithm section see the algorithm's documentation: @n
    @ref OpenMS::PrecursorIonSelection @n

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPrecursorIonSelector :
  public TOPPBase
{
public:
  TOPPPrecursorIonSelector() :
    TOPPBase("PrecursorIonSelector", "PrecursorIonSelector")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<input file>", "", "Input feature map file (featureXML)");
    setValidFormats_("in", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out", "<output file>", "", "modified feature map", false);
    setValidFormats_("out", ListUtils::create<String>("featureXML"));

    registerOutputFile_("next_feat", "<output file>", "", "feature map (featureXML) file with the selected precursors", false);
    setValidFormats_("next_feat", ListUtils::create<String>("featureXML"));

    registerInputFile_("ids", "<id file>", "", "file containing results of identification");
    setValidFormats_("ids", ListUtils::create<String>("idXML"));
    
    registerIntOption_("num_precursors", "<Int>", 1, "number of precursors to be selected", false);
    registerInputFile_("raw_data", "<file>", "", "Input profile data.", false);
    setValidFormats_("raw_data", ListUtils::create<String>("mzML"));
    registerFlag_("load_preprocessing", "The preprocessed db is loaded from file, not calculated.");
    registerFlag_("store_preprocessing", "The preprocessed db is stored.");
    registerFlag_("simulation", "Simulate the whole LC-MS/MS run.");
    registerOutputFile_("sim_results", "<output file>", "", "File containing the results of the simulation run", false);
    setValidFormats_("sim_results", ListUtils::create<String>("txt"));
    
    registerInputFile_("db_path", "<db-file>", "", "db file", false);
    setValidFormats_("db_path", ListUtils::create<String>("fasta"));

    registerInputFile_("rt_model", "<rt-model-file>", "", "SVM Model for RTPredict", false);
    setValidFormats_("rt_model", ListUtils::create<String>("txt"));
        
    registerInputFile_("dt_model", "<dt-model-file>", "", "SVM Model for PTPredict", false);
    setValidFormats_("dt_model", ListUtils::create<String>("txt"));
    
    registerStringOption_("solver", "<solver-type>", "GLPK", "LP solver type", false, true);
    setValidStrings_("solver", ListUtils::create<String>("GLPK,COINOR"));
    registerStringList_("fixed_modifications", "<mods>", StringList(), "the modifications i.e. Carboxymethyl (C)", false);
    addEmptyLine_();
    registerSubsection_("algorithm", "Settings for the compound list creation and rescoring.");
  }

  Param getSubsectionDefaults_(const String & /* section*/) const override
  {
    Param param = PrecursorIonSelection().getDefaults();
    //    param.insert("",PrecursorIonSelection().getDefaults().copy(""));
    return param;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String in(getStringOption_("in"));
    String out(getStringOption_("out"));
    String raw_data(getStringOption_("raw_data"));
    String next_prec = getStringOption_("next_feat");
    String ids = getStringOption_("ids");
    String db_path = getStringOption_("db_path");
    UInt prec_num = getIntOption_("num_precursors");
    bool simulation = getFlag_("simulation");
    String sim_results = getStringOption_("sim_results");
    bool load_preprocessing = getFlag_("load_preprocessing");
    bool store_preprocessing = getFlag_("store_preprocessing");
    String rt_model = getStringOption_("rt_model");
    String dt_model = getStringOption_("dt_model");
    String solver(getStringOption_("solver"));
    StringList fixed_mods = getStringList_("fixed_modifications");
    //-------------------------------------------------------------
    // init pis preprocessing
    //-------------------------------------------------------------
    Param pisp_param = getParam_().copy("algorithm:Preprocessing:", true);
    pisp_param.remove("type");
    pisp_param.remove("min_pep_ids");
    pisp_param.remove("max_iteration");
    writeDebug_("Parameters passed to PrecursorIonSelectionPreprocessing", pisp_param, 3);
    PrecursorIonSelectionPreprocessing pisp;
    //    pisp.setLogType(log_type_);
    pisp.setParameters(pisp_param);
    pisp.setFixedModifications(fixed_mods);
    if (load_preprocessing)
    {
      pisp.loadPreprocessing();
    }
    else if (db_path == "")
    {
      writeLog_("No database file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    else if (rt_model == "" ||  dt_model == "")
    {
      pisp.dbPreprocessing(db_path, store_preprocessing);
    }
    else
    {
      pisp.dbPreprocessing(db_path, rt_model, dt_model, store_preprocessing);
    }

    PeakMap exp;
    if (raw_data != "") MzMLFile().load(raw_data, exp);

    //-------------------------------------------------------------
    // init pis
    //-------------------------------------------------------------
    Param pis_param = getParam_().copy("algorithm:", true);
    pis_param.removeAll("preprocessing");
    writeDebug_("Parameters passed to PrecursorIonSelection", pis_param, 3);
    PrecursorIonSelection pis;
    //    pis.setLogType(log_type_);
    pis.setParameters(pis_param);
#if COINOR_SOLVER == 1
    if (solver == "GLPK")
    {
      pis.setLPSolver(LPWrapper::SOLVER_GLPK);
    }
    else pis.setLPSolver(LPWrapper::SOLVER_COINOR);
#endif

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    FeatureMap f_map;
    FeatureXMLFile f_file;
    f_file.load(in, f_map);

    std::vector<PeptideIdentification> pep_ids;
    std::vector<ProteinIdentification> prot_ids;
    String document_id;
    IdXMLFile idxml_file;
    idxml_file.load(ids, prot_ids, pep_ids, document_id);

    //-------------------------------------------------------------
    // preprocessing, rescoring
    //-------------------------------------------------------------

    if (simulation)
    {
      pis.simulateRun(f_map, pep_ids, prot_ids, pisp, sim_results, exp, "");
    }
    else
    {

      pis.rescore(f_map, pep_ids, prot_ids, pisp);  // todo: add "rescoring" for LP selection
      FeatureMap new_precursors;
      pis.getNextPrecursors(f_map, new_precursors, prec_num);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      if (next_prec != "") f_file.store(next_prec, new_precursors);
    }

    if (out != "") f_file.store(out, f_map);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPPrecursorIonSelector tool;
  return tool.main(argc, argv);
}

/// @endcond
