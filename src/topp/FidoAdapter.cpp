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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QDir>

#include <fstream>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_FidoAdapter FidoAdapter

    @brief Runs the protein inference engine Fido.

<CENTER>
    <table>
        <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FidoAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_PeptideIndexer\n(with @p annotate_proteins option) </td>
            <td VALIGN="middle" ALIGN="center" ROWSPAN=2> @ref TOPP_ProteinQuantifier\n(via @p protein_groups parameter) </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability\n(with @p prob_correct option) </td>
        </tr>
    </table>
</CENTER>

    This tool wraps the protein inference algorithm Fido (http://noble.gs.washington.edu/proj/fido/). Fido uses a Bayesian probabilistic model to group and score proteins based on peptide-spectrum matches. It was published in:

    Serang <em>et al.</em>: <a href="http://pubs.acs.org/doi/abs/10.1021/pr100594k">Efficient marginalization to compute protein posterior probabilities from shotgun mass spectrometry data</a> (J. Proteome Res., 2010).

    <b>Input format:</b>

    Care has to be taken to provide suitable input data for this adapter. In the peptide/protein identification results (e.g. coming from a database search engine), the proteins have to be annotated with target/decoy meta data. To achieve this, run @ref TOPP_PeptideIndexer with the @p annotate_proteins option switched on.@n
    In addition, the scores for peptide hits in the input data have to be posterior probabilities - as produced e.g. by PeptideProphet in the TPP or by @ref TOPP_IDPosteriorErrorProbability (with the @p prob_correct option switched on) in OpenMS. Inputs from @ref TOPP_IDPosteriorErrorProbability (without @p prob_correct) or from @ref TOPP_ConsensusID are treated as special cases: Their posterior error probabilities (lower is better) are converted to posterior probabilities (higher is better) for processing.

    <b>Output format:</b>

    The output of this tool is an augmented version of the input: The protein groups and accompanying posterior probabilities inferred by Fido are stored as "indistinguishable protein groups", attached to the (first) protein identification run of the input data.@n
    The result can be passed to @ref TOPP_ProteinQuantifier via its @p protein_groups parameter, to have the protein grouping taken into account during quantification.


    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_FidoAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_FidoAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFidoAdapter :
  public TOPPBase
{
public:
  TOPPFidoAdapter() :
    TOPPBase("FidoAdapter", "Runs the protein inference engine Fido.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input: identification results");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output: identification results "
                        "with scored/grouped proteins");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerStringOption_("exe", "<file>", 
#if defined(OPENMS_WINDOWSPLATFORM)
                          "FidoChooseParameters.exe",
#else
                          "FidoChooseParameters",
#endif
                          "Executable for Fido with parameter estimation", 
                          false);

    registerFlag_("no_cleanup", "Omit clean-up of peptide sequences (removal of non-letter characters, replacement of I with L)");
    registerFlag_("all_PSMs", "Consider all PSMs of each peptide, instead of only the best one");
    registerFlag_("group_level", "Perform inference on protein group level (instead of individual protein level). This will lead to higher probabilities for (bigger) protein groups.");
    registerStringOption_("accuracy", "<choice>", "", "Accuracy level of start parameters. There is a trade-off between accuracy and runtime. Empty uses the default ('best').", false, true);
    setValidStrings_("accuracy", ListUtils::create<String>(",best,relaxed,sloppy"));
    registerIntOption_("log2_states", "<number>", 0, "Binary logarithm of the max. number of connected states in a subgraph. For a value N, subgraphs that are bigger than 2^N will be split up, sacrificing accuracy for runtime. '0' uses the default (18).", false);
    setMinInt_("log2_states", 0);
    registerIntOption_("log2_states_precalc", "<number>", 0, "Like 'log2_states', but allows to set a separate limit for the precalculation", false, true);
    setMinInt_("log2_states_precalc", 0);
  }


  ExitCodes main_(int, const char**)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String exe = getStringOption_("exe");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;

    LOG_INFO << "Reading input data..." << endl;
    IdXMLFile().load(in, proteins, peptides);
    if (proteins.empty() || peptides.empty())
    {
      LOG_ERROR << "Error: Input file '" << in 
                << "' should contain both protein and peptide data." << endl;
      return INPUT_FILE_EMPTY;
    }
    if (proteins.size() > 1)
    {
      LOG_WARN << "Warning: Input contains more than one protein "
               << "identification run. Protein inference results will be "
               << "written to the first run only." << endl;
    }

    LOG_INFO << "Generating temporary files for Fido..." << endl;
    String temp_directory = 
      QDir::toNativeSeparators((File::getTempDirectory() + "/" + 
                                File::getUniqueName() + "/").toQString());
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }

    String fido_input_graph = temp_directory + "fido_input_graph.txt";
    String fido_input_proteins = temp_directory + "fido_input_proteins.txt";

    // write PSM graph:
    ofstream graph_out(fido_input_graph.c_str());
    bool warned_once = false;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty()) continue;
      pep_it->sort();
      const PeptideHit& hit = pep_it->getHits()[0];
      if (hit.getSequence().empty() || 
          hit.getProteinAccessions().empty()) continue;
      double score = hit.getScore();

      String error_reason;
      if (!pep_it->isHigherScoreBetter())
      {
        // workaround for important TOPP tools:
        String score_type = pep_it->getScoreType();
        score_type.toLower();
        if ((score_type == "posterior error probability") || 
            score_type.hasPrefix("consensus_"))
        {
          if (!warned_once)
          {
            LOG_WARN << "Warning: Scores of peptide hits seem to be posterior "
              "error probabilities. Converting to (positive) posterior "
              "probabilities." << endl;
            warned_once = true;
          }
          score = 1.0 - score;
        }
        else error_reason = "lower scores are better";
      }
      else if (score < 0.0)
      {
        error_reason = "score < 0";
      }
      else if (score > 1.0)
      {
        error_reason = "score > 1";
      }
      if (!error_reason.empty())
      {
        String msg = "Error: Unsuitable score type for peptide-spectrum "
          "matches detected (problem: " + error_reason + ").\nFido requires "
          "probabilities as scores, e.g. as produced by "
          "IDPosteriorErrorProbability with the 'prob_correct' option.";
        LOG_ERROR << msg << endl;
        return INCOMPATIBLE_INPUT_DATA;
      }

      graph_out << "e " << hit.getSequence() << endl;
      for (vector<String>::const_iterator acc_it = 
             hit.getProteinAccessions().begin(); acc_it != 
             hit.getProteinAccessions().end(); ++acc_it)
      {
        graph_out << "r " << *acc_it << endl;
      }
      graph_out << "p " << score << endl;
    }
    graph_out.close();

    // gather protein target/decoy data:
    set<String> targets, decoys;
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      for (vector<ProteinHit>::iterator hit_it = prot_it->getHits().begin();
           hit_it != prot_it->getHits().end(); ++hit_it)
      {
        String target_decoy = hit_it->getMetaValue("target_decoy").toString();
        if (target_decoy == "target")
        {
          targets.insert(hit_it->getAccession());
        }
        else if (target_decoy == "decoy")
        {
          decoys.insert(hit_it->getAccession());
        }
        else
        {
          String msg = "Error: All protein hits must be annotated with target/"
            "decoy meta data. Run PeptideIndexer with the 'annotate_proteins' "
            "option to accomplish this.";
          LOG_ERROR << msg << endl;
          return INCOMPATIBLE_INPUT_DATA;
        }
      }
    }

    // write target/decoy protein sets:
    ofstream proteins_out(fido_input_proteins.c_str());
    proteins_out << "{ ";
    for (set<String>::iterator it = targets.begin(); it != targets.end(); ++it)
    {
      if (it != targets.begin()) proteins_out << " , ";
      proteins_out << *it;
    }
    proteins_out << " }\n{ ";
    for (set<String>::iterator it = decoys.begin(); it != decoys.end(); ++it)
    {
      if (it != decoys.begin()) proteins_out << " , ";
      proteins_out << *it;
    }
    proteins_out << " }" << endl;
    proteins_out.close();

    LOG_INFO << "Running Fido..." << endl;
    // Fido parameters:
    QStringList inputs;
    if (getFlag_("no_cleanup")) inputs << "-p";
    if (getFlag_("all_PSMs")) inputs << "-a";
    if (getFlag_("group_level")) inputs << "-g";
    String accuracy = getStringOption_("accuracy");
    if (!accuracy.empty())
    {
      if (accuracy == "best") inputs << "-c 1";
      else if (accuracy == "relaxed") inputs << "-c 2";
      else if (accuracy == "sloppy") inputs << "-c 3";
    }
    inputs << fido_input_graph.toQString() << fido_input_proteins.toQString();
    Int log2_states = getIntOption_("log2_states");
    Int log2_states_precalc = getIntOption_("log2_states_precalc");
    if (log2_states_precalc)
    {
      if (!log2_states) log2_states = 18; // actual default value
      inputs << QString::number(log2_states_precalc);
    }
    if (log2_states) inputs << QString::number(log2_states);

    // run program and read output:
    QProcess fido;
    fido.start(exe.toQString(), inputs);

    ExitCodes exit_code = EXECUTION_OK;
    if (!fido.waitForFinished(-1))
    {
      String msg = "Fatal error running Fido (command: '" + exe + " \"" + 
        String(inputs.join("\" \"")) + "\"').\nDoes the Fido executable exist?";
      LOG_ERROR << msg << endl;
      exit_code = EXTERNAL_PROGRAM_ERROR;
    }
    else // success!
    {
      LOG_INFO << "Parsing Fido results and writing output..." << endl;
      String output = QString(fido.readAllStandardOutput());
      if (debug_level_ > 1)
      {
        String fido_output = temp_directory + "fido_output.txt";
        ofstream results(fido_output.c_str());
        results << output;
        results.close();
      }

      vector<String> lines;
      output.split("\n", lines);

      Int protein_counter = 0;
      vector<ProteinIdentification::ProteinGroup> groups;
      for (vector<String>::iterator line_it = lines.begin(); 
           line_it != lines.end(); ++line_it)
      {
        // format of a line (example):
        // 0.6788 { SW:TRP6_HUMAN , GP:AJ271067_1 , GP:AJ271068_1 }
        istringstream line(*line_it);
        ProteinIdentification::ProteinGroup group;
        line >> group.probability;
        // parse accessions (won't work if accessions can contain spaces!):
        String accession;
        line >> accession;
        while (line)
        {
          if (accession.size() > 1) // skip braces and commas
          {
            group.accessions.push_back(accession);
          }
          line >> accession;
        }
        if (!group.accessions.empty())
        {
          protein_counter += group.accessions.size();
          groups.push_back(group);
        }
      }
      proteins[0].getIndistinguishableProteins() = groups;
      LOG_INFO << "Inferred " << protein_counter << " proteins in "
               << groups.size() << " groups." << endl;

      // write output:
      IdXMLFile().store(out, proteins, peptides);
    }

    // clean up temporary files
    if (debug_level_ > 1)
    {
      LOG_DEBUG << "Keeping temporary files at '" << temp_directory
                << "'. Set debug level to 0 or 1 to remove them." << endl;
    }
    else
    {
      LOG_INFO << "Removing temporary files..." << endl;
      File::removeDirRecursively(temp_directory);
      if (debug_level_ == 1)
      {
        String msg = "Set debug level to 2 or higher to keep temporary files "
          "at '" + temp_directory + "'.";
        LOG_DEBUG << msg << endl;
      }
    }

    return exit_code;
  }

};


int main(int argc, const char** argv)
{
  TOPPFidoAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
