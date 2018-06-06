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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Xiao Liang, Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>

#include <boost/bimap.hpp>

#include <QtCore/QProcess>
#include <QDir>

#include <fstream>
#include <sstream>
#include <queue>

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
       <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ FidoAdapter \f$ \longrightarrow \f$</td>
       <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
     </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_PeptideIndexer </td>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=3> @ref TOPP_ProteinQuantifier\n(via @p protein_groups parameter) </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability\n(with @p prob_correct option) </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref UTILS_IDScoreSwitcher </td>
    </tr>
   </table>
 </CENTER>

 This tool wraps the protein inference algorithm Fido (http://noble.gs.washington.edu/proj/fido/). Fido uses a Bayesian probabilistic model to group and score proteins based on peptide-spectrum matches. It was published in:

 Serang <em>et al.</em>: <a href="http://pubs.acs.org/doi/abs/10.1021/pr100594k">Efficient marginalization to compute protein posterior probabilities from shotgun mass spectrometry data</a> (J. Proteome Res., 2010).

 By default, this adapter runs the Fido variant with parameter estimation (@p FidoChooseParameters), as recommended by the authors of Fido. However, it is also possible to run "pure" Fido by setting the @p prob:protein, @p prob:peptide and @p prob:spurious parameters, if appropriate values are known (e.g. from a previous Fido run). Other parameters, except for @p log2_states, are not applicable in this case.

 Depending on the @p separate_runs setting, data from input files containing multiple protein identification runs (e.g. several replicates or different search engines) will be merged (default) or annotated separately.

 <b>Input format:</b>

 Care has to be taken to provide suitable input data for this adapter. In the peptide/protein identification results (e.g. coming from a database search engine), the proteins have to be annotated with target/decoy meta data. To achieve this, run @ref TOPP_PeptideIndexer.@n
 In addition, the scores for peptide hits in the input data have to be posterior probabilities - as produced e.g. by PeptideProphet in the TPP or by @ref TOPP_IDPosteriorErrorProbability (with the @p prob_correct option switched on) in OpenMS. If scores are found to be posterior error probabilities (PEPs, lower is better), they are converted to posterior probabilities (higher is better) using "1 - PEP".@n
 If the posterior (error) probabilities are stored in user parameters ("UserParam") in the idXML instead of in the score fields, @ref UTILS_IDScoreSwitcher can be used to rewrite the scores. (This may be the case e.g. if @ref TOPP_FalseDiscoveryRate and @ref TOPP_IDFilter were applied for FDR filtering prior to protein inference.)

 <b>Output format:</b>

 The output of this tool is an augmented version of the input: The protein groups and accompanying posterior probabilities inferred by Fido are stored as "indistinguishable protein groups", attached to the protein identification run(s) of the input data. Also attached are meta values recording the Fido parameters (@p Fido_prob_protein, @p Fido_prob_peptide, @p Fido_prob_spurious).@n
 The result can be passed to @ref TOPP_ProteinQuantifier via its @p protein_groups parameter, to have the protein grouping taken into account during quantification.@n
 Note that if the input contains multiple identification runs and @p separate_runs is @e not set (the default), the identification data from all runs will be pooled for the Fido analysis and the result will only contain one (merged) identification run. This is the desired behavior if the protein grouping should be used by ProteinQuantifier.
 When the @p greedy_group_resolution flag is set, "peptide to indistinguishable proteins" mappings will be unique in the output and the actual resolved groups are added as "protein groups", attached to the protein identification run(s) of the input data (in addition to the "indistinguishable protein groups").

 @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

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

  typedef boost::bimap<String, String> StringBimap;

  StringBimap sanitized_accessions_; // protein accessions

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input: identification results");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output: identification results with scored/grouped proteins");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("fido_executable", "<path>", "Fido", "Path to the Fido executable to use; may be empty if the executable is globally available.", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("fidocp_executable", "<path>", "FidoChooseParameters", "Path to the FidoChooseParameters executable to use; may be empty if the executable is globally available.", true, false, ListUtils::create<String>("skipexists"));
    registerFlag_("separate_runs", "Process multiple protein identification runs in the input separately, don't merge them. Merging results in loss of descriptive information of the single protein identification runs.");
    registerFlag_("keep_zero_group", "Keep the group of proteins with estimated probability of zero, which is otherwise removed (it may be very large)", true);
    registerFlag_("greedy_group_resolution", "Post-process Fido output with greedy resolution of shared peptides based on the protein probabilities. Also adds the resolved ambiguity groups to output.");
    registerFlag_("no_cleanup", "Omit clean-up of peptide sequences (removal of non-letter characters, replacement of I with L)");
    registerFlag_("all_PSMs", "Consider all PSMs of each peptide, instead of only the best one");
    registerFlag_("group_level", "Perform inference on protein group level (instead of individual protein level). This will lead to higher probabilities for (bigger) protein groups.");
    registerStringOption_("accuracy", "<choice>", "", "Accuracy level of start parameters. There is a trade-off between accuracy and runtime. Empty uses the default ('best').", false, true);
    setValidStrings_("accuracy", ListUtils::create<String>(",best,relaxed,sloppy"));
    registerIntOption_("log2_states", "<number>", 0, "Binary logarithm of the max. number of connected states in a subgraph. For a value N, subgraphs that are bigger than 2^N will be split up, sacrificing accuracy for runtime. '0' uses the default (18).", false);
    setMinInt_("log2_states", 0);
    registerIntOption_("log2_states_precalc", "<number>", 0, "Like 'log2_states', but allows to set a separate limit for the precalculation", false, true);
    setMinInt_("log2_states_precalc", 0);
    registerTOPPSubsection_("prob", "Probability values for running Fido directly, i.e. without parameter estimation (in which case other settings, except 'log2_states', are ignored)");
    registerDoubleOption_("prob:protein", "<value>", 0.0, "Protein prior probability ('gamma' parameter)", false);
    setMinFloat_("prob:protein", 0.0);
    registerDoubleOption_("prob:peptide", "<value>", 0.0, "Peptide emission probability ('alpha' parameter)", false);
    setMinFloat_("prob:peptide", 0.0);
    registerDoubleOption_("prob:spurious", "<value>", 0.0, "Spurious peptide identification probability ('beta' parameter)", false);
    setMinFloat_("prob:spurious", 0.0);
  }

  // write a PSM graph file for Fido based on the given peptide identifications;
  // optionally only use IDs with given identifier (filter by protein ID run):
  void writePSMGraph_(vector<PeptideIdentification>& peptides,
                      const String& out_path, const String& identifier = "")
  {
    ofstream graph_out(out_path.c_str());
    bool warned_once = false;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if ((!identifier.empty() && (pep_it->getIdentifier() != identifier)) ||
          (pep_it->getHits().empty()))
      {
        continue;
      }
      pep_it->sort();
      const PeptideHit& hit = pep_it->getHits()[0];
      if (hit.getSequence().empty() || hit.extractProteinAccessionsSet().empty())
      {
        continue;
      }

      double score = hit.getScore();
      String score_type = pep_it->getScoreType();
      bool higher_better = pep_it->isHigherScoreBetter();

      // workaround for posterior error probabilities (PEPs):
      if (score_type.hasSuffix("_score"))
      {
        score_type.chop(6);
      }
      score_type.toLower().remove(' ').remove('_').remove('-');
      if (score_type == "posteriorerrorprobability")
      {
        if (!warned_once)
        {
          LOG_WARN << "Warning: Scores of peptide hits seem to be posterior "
            "error probabilities. Converting to (positive) posterior "
            "probabilities." << endl;
          warned_once = true;
        }
        score = 1.0 - score;
        higher_better = true;
      }

      String error_reason;
      if (!higher_better)
      {
        error_reason = "lower scores are better";
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
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      }

      // Remove modifications before writing to input graph file
      graph_out << "e " << hit.getSequence().toUnmodifiedString() << endl;
      const set<String>& accessions = hit.extractProteinAccessionsSet();
      for (set<String>::const_iterator acc_it = accessions.begin();
           acc_it != accessions.end(); ++acc_it)
      {
        if (acc_it->empty()) continue;
        graph_out << "r " << sanitized_accessions_.left.find(*acc_it)->second
                  << endl;
      }
      graph_out << "p " << score << endl;
    }
    graph_out.close();
  }


  // write the list of target and decoy proteins for Fido:
  void writeProteinLists_(const ProteinIdentification& protein,
                          const String& out_path)
  {
    // gather protein target/decoy data:
    set<String> targets, decoys;
    for (vector<ProteinHit>::const_iterator hit_it = protein.getHits().begin();
         hit_it != protein.getHits().end(); ++hit_it)
    {
      String target_decoy = hit_it->getMetaValue("target_decoy").toString();
      String accession = hit_it->getAccession();
      String sanitized = sanitized_accessions_.left.find(accession)->second;
      if (target_decoy == "target")
      {
        targets.insert(sanitized);
      }
      else if (target_decoy == "decoy")
      {
        decoys.insert(sanitized);
      }
      else
      {
        String msg = "Error: All protein hits must be annotated with target/"
          "decoy meta data. Run PeptideIndexer with the 'annotate_proteins' "
          "option to accomplish this.";
        LOG_ERROR << msg << endl;
        throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    if (targets.empty())
    {
      String msg = "Error: No target proteins found. Fido needs both targets "
        "and decoys.";
      LOG_ERROR << msg << endl;
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    if (decoys.empty())
    {
      String msg = "Error: No decoy proteins found. Fido needs both targets "
        "and decoys.";
      LOG_ERROR << msg << endl;
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
    }

    // write sets to file:
    ofstream proteins_out(out_path.c_str());
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
  }


  // run Fido(ChooseParameters) and read output:
  bool runFido_(ProteinIdentification& protein,
                vector<PeptideIdentification>& peptides, bool choose_params,
                const String& exe, QStringList& fido_params,
                double& prob_protein, double& prob_peptide,
                double& prob_spurious, const String& temp_dir,
                bool keep_zero_group = false, bool greedy_flag = false,
                Size counter = 0)
  {
    // create a copy of the params so the templates can be overwritten with
    // different values:
    QStringList current_fido_params = QStringList(fido_params);

    LOG_INFO << "Generating temporary files for Fido..." << endl;
    String num = counter ? "." + String(counter) : "";
    String input_graph = temp_dir + "fido_input_graph" + num + ".txt";
    current_fido_params.replaceInStrings("INPUT_GRAPH",
                                         input_graph.toQString());

    writePSMGraph_(peptides, input_graph, protein.getIdentifier());
    if (choose_params)
    {
      String input_proteins = temp_dir + "fido_input_proteins" + num + ".txt";
      current_fido_params.replaceInStrings("INPUT_PROTEINS",
                                           input_proteins.toQString());
      writeProteinLists_(protein, input_proteins);
      LOG_INFO << "Running Fido with parameter estimation..." << endl;
    }
    else LOG_INFO << "Running Fido with fixed parameters..." << endl;

    QProcess fido;
    fido.start(exe.toQString(), current_fido_params);

    if (!fido.waitForFinished(-1))
    {
      String cmd = exe + " \"" + String(current_fido_params.join("\" \"")) +
        "\"";
      LOG_ERROR << "Fatal error running Fido (command: '" + cmd + "').\n"
                << "Does the Fido executable exist?" << endl;
      return false;
    }

    // success! parse output:
    vector<String> lines;
    if (choose_params) // get relevant parts of parameter search output
    {
      String params_output = QString(fido.readAllStandardError());
      LOG_INFO << "Fido parameter search:" << endl;
      if (debug_level_ > 1)
      {
        String output_status = temp_dir + "fido_status" + num + ".txt";
        ofstream status(output_status.c_str());
        status << params_output;
        status.close();
      }
      params_output.split("\n", lines);
      vector<String>::iterator pos = remove(lines.begin(), lines.end(), "");
      lines.erase(pos, lines.end());
      if (!lines.empty())
      {
        if (lines[0].hasPrefix("caught an exception"))
        {
          LOG_ERROR << "Error running Fido: '" + lines[0] + "'" << endl;
          return false;
        }
        if (lines[0].hasPrefix("Warning:")) LOG_WARN << lines[0] << endl;
        if (lines.back().hasPrefix("Using best gamma, alpha, beta ="))
        {
          LOG_INFO << lines.back() << endl;
          stringstream ss;
          ss << lines.back().suffix('=');
          ss >> prob_protein >> prob_peptide >> prob_spurious;
        }
      }
    }

    LOG_INFO << "Parsing Fido results and writing output..." << endl;
    String output = QString(fido.readAllStandardOutput());
    if (debug_level_ > 1)
    {
      String output_result = temp_dir + "fido_output" + num + ".txt";
      ofstream results(output_result.c_str());
      results << output;
      results.close();
    }
    output.split("\n", lines);

    // edit ProteinIdentification with Fido information and scores but save old
    // information and scores in MetaInfo:
    const String old_scoretype = protein.getScoreType();
    protein.setScoreType("Posterior Probability");
    protein.setHigherScoreBetter(true);
    protein.setDateTime(DateTime::now());

    Size protein_counter = 0, zero_proteins = 0;

    // add protein groups to ProteinIdentification:
    vector<ProteinIdentification::ProteinGroup> groups;
    for (vector<String>::iterator line_it = lines.begin();
         line_it != lines.end(); ++line_it)
    {
      // format of a line (example):
      // 0.6788 { SW:TRP6_HUMAN , GP:AJ271067_1 , GP:AJ271068_1 }
      istringstream line(*line_it);
      ProteinIdentification::ProteinGroup group;
      line >> group.probability;

      // parse accessions:
      String accession;
      while (line)
      {
        line >> accession;
        if (accession.size() > 1) // skip braces and commas
        {
          if (group.probability == 0.0)
          {
            ++zero_proteins;
            if (!keep_zero_group) continue;
          }

          // de-sanitize:
          accession = sanitized_accessions_.right.find(accession)->second;

          // look up corresponding hits and update scores...
          std::vector<ProteinHit>::iterator hit = protein.findHit(accession);

          // ...while saving possible old ones (if there are any):
          if (!old_scoretype.empty())
          {
            hit->setMetaValue(old_scoretype + "_score", hit->getScore());
          }
          else if (hit->metaValueExists("old_score_type"))
          {
            // save old score type:
            String old_scoretype_meta = hit->getMetaValue("old_score_type");
            if (!old_scoretype_meta.empty())
            {
              hit->setMetaValue(old_scoretype_meta + "_score", hit->getScore());
            }
            hit->removeMetaValue("old_score_type");
          }

          hit->setScore(group.probability);
          group.accessions.push_back(accession);
        }
      }
      if (!group.accessions.empty())
      {
        protein_counter += group.accessions.size();
        sort(group.accessions.begin(), group.accessions.end());
        groups.push_back(group);
      }
    }

    // Sort groups by probability and add the finally used Fido params
    // as meta values
    sort(groups.begin(), groups.end());
    protein.getIndistinguishableProteins() = groups;
    protein.setMetaValue("Fido_prob_protein", prob_protein);
    protein.setMetaValue("Fido_prob_peptide", prob_peptide);
    protein.setMetaValue("Fido_prob_spurious", prob_spurious);
    LOG_INFO << "Inferred " << protein_counter << " proteins in "
             << groups.size() << " groups ("
             << ((keep_zero_group && zero_proteins) ? "including " : "")
             << zero_proteins << " proteins with probability zero"
             << ((keep_zero_group || !zero_proteins) ? ")." : " not included).")
             << endl;

    // Do post-processing on groups if specified
    if (greedy_flag)
    {
      LOG_INFO << "Resolving ambiguity groups greedily on Fido output..."
               << endl;
      PeptideProteinResolution graph = PeptideProteinResolution();
      graph.buildGraph(protein, peptides);
      graph.resolveGraph(protein, peptides);
    }
    return true;
  }


  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String fido_executable = getStringOption_("fido_executable");
    String fidocp_executable = getStringOption_("fidocp_executable");
    bool separate_runs = getFlag_("separate_runs");
    bool keep_zero_group = getFlag_("keep_zero_group");
    bool greedy_flag = getFlag_("greedy_group_resolution");
    double prob_protein = getDoubleOption_("prob:protein");
    double prob_peptide = getDoubleOption_("prob:peptide");
    double prob_spurious = getDoubleOption_("prob:spurious");
    bool choose_params = ((prob_protein == 0.0) && (prob_peptide == 0.0) &&
                          (prob_spurious == 0.0)); // use FidoChooseParameters?

    if (fido_executable.empty()) // expect executables in PATH
    {
      fido_executable = "Fido";
    }

    if (fidocp_executable.empty()) // expect executables in PATH
    {
      fidocp_executable = "FidoChooseParameters";
    }
    String executable = choose_params ? fidocp_executable : fido_executable;

    // input data:
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

    // remove protein hits that shouldn't be there:
    IDFilter::removeUnreferencedProteins(proteins, peptides);

    // sanitize protein accessions:
    set<String> accessions;
    for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
         prot_it != proteins.end(); ++prot_it)
    {
      for (vector<ProteinHit>::iterator hit_it = prot_it->getHits().begin();
           hit_it != prot_it->getHits().end(); ++hit_it)
      {
        accessions.insert(hit_it->getAccession());
      }
    }
    Size acc_counter = 1;
    for (set<String>::iterator it = accessions.begin(); it != accessions.end();
         ++it, ++acc_counter)
    {
      // take valid prefix (= accession) and add number to ensure uniqueness:
      Size pos = it->find_first_of(" \t,{}");
      // using "it->prefix(pos)" doesn't work if none of the chars was found:
      String sanitized = it->substr(0, pos) + "_" + String(acc_counter);
      sanitized_accessions_.insert(StringBimap::value_type(*it, sanitized));
    }

    // create temporary directory:
    String temp_dir = File::getTempDirectory() + "/" + File::getUniqueName() +
      "/";
    temp_dir = QDir::toNativeSeparators(temp_dir.toQString());
    {
      QDir d;
      d.mkpath(temp_dir.toQString());
    }

    // Fido parameters (use placeholders for paths - replace them later):
    QStringList fido_params;
    Int log2_states = getIntOption_("log2_states");
    if (choose_params)
    {
      if (getFlag_("no_cleanup")) fido_params << "-p";
      if (getFlag_("all_PSMs")) fido_params << "-a";
      if (getFlag_("group_level")) fido_params << "-g";
      String accuracy = getStringOption_("accuracy");
      if (!accuracy.empty())
      {
        if (accuracy == "best") fido_params << "-c 1";
        else if (accuracy == "relaxed") fido_params << "-c 2";
        else if (accuracy == "sloppy") fido_params << "-c 3";
      }
      fido_params << "INPUT_GRAPH" << "INPUT_PROTEINS";
      Int log2_states_precalc = getIntOption_("log2_states_precalc");
      if (log2_states_precalc)
      {
        if (!log2_states) log2_states = 18; // actual default value
        fido_params << QString::number(log2_states_precalc);
      }
    }
    else // run Fido only
    {
      fido_params << "INPUT_GRAPH" << QString::number(prob_protein)
                  << QString::number(prob_peptide)
                  << QString::number(prob_spurious);
    }
    if (log2_states) fido_params << QString::number(log2_states);


    // actually run Fido now and process its output:
    bool fido_success = false;

    if (separate_runs || proteins.size() == 1)
    {
      // treat multiple protein ID runs separately or process the only run:
      Size counter = 1;
      for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
           prot_it != proteins.end(); ++prot_it, ++counter)
      {
        LOG_INFO << "Protein identification run " << counter << ":" << endl;
        fido_success = runFido_(*prot_it, peptides, choose_params, executable,
                                fido_params, prob_protein, prob_peptide,
                                prob_spurious, temp_dir, keep_zero_group,
                                greedy_flag, counter);
      }
    }
    else // merge multiple protein ID runs
    {
      ProteinIdentification all_proteins; // one ID run to merge all hits
      // set search engine to Fido since they might disagree for different runs:
      all_proteins.setSearchEngine("Fido");

      // make sure identifiers match (otherwise "IdXMLFile::store" complains):
      all_proteins.setIdentifier("");
      for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        pep_it->setIdentifier("");
      }

      // for every protein (accession), save the first occurrence:
      map<String, ProteinHit*> hit_map; // protein hits by accession
      for (vector<ProteinIdentification>::reverse_iterator prot_it =
           proteins.rbegin(); prot_it != proteins.rend(); ++prot_it)
      {
        for (vector<ProteinHit>::reverse_iterator hit_it =
               prot_it->getHits().rbegin(); hit_it !=
               prot_it->getHits().rend(); ++hit_it)
        {
          // save original score type:
          // @TODO: does this make sense if we have potentially different score
          // types from different search engines?
          hit_it->setMetaValue("old_score_type", prot_it->getScoreType());
          hit_map[hit_it->getAccession()] = &(*hit_it);
        }
      }
      all_proteins.getHits().reserve(hit_map.size());
      for (map<String, ProteinHit*>::iterator map_it = hit_map.begin();
           map_it != hit_map.end(); ++map_it)
      {
        all_proteins.insertHit(*(map_it->second));
      }

      fido_success = runFido_(all_proteins, peptides, choose_params, executable,
                              fido_params, prob_protein, prob_peptide,
                              prob_spurious, temp_dir, keep_zero_group,
                              greedy_flag);

      // replace proteins with merged variant:
      proteins.clear();
      proteins.push_back(all_proteins);
    }

    // write output:
    IdXMLFile().store(out, proteins, peptides);

    // clean up temporary files:
    if (debug_level_ > 1)
    {
      LOG_DEBUG << "Keeping temporary files at '" << temp_dir
                << "'. Set debug level to 0 or 1 to remove them." << endl;
    }
    else
    {
      LOG_INFO << "Removing temporary files..." << endl;
      File::removeDirRecursively(temp_dir);
      if (debug_level_ == 1)
      {
        String msg = "Set debug level to 2 or higher to keep temporary files "
          "at '" + temp_dir + "'.";
        LOG_DEBUG << msg << endl;
      }
    }

    return fido_success ? EXECUTION_OK : EXTERNAL_PROGRAM_ERROR;
  }
};

int main(int argc, const char** argv)
{
  TOPPFidoAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
