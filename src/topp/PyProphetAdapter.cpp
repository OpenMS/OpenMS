// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QStringList>

#include <cctype>
#include <unordered_set>

/**
@page TOPP_PyProphetAdapter PyProphetAdapter
@brief Adapter to run PyProphet scoring, optional peptide/protein/gene inference, peptidoform inference (IPF), and exports on OSW files.
\see http://openswath.org/

PyProphetAdapter orchestrates:\n
 - (optional) \em merge of multiple OSW inputs (\c pyprophet merge osw)\n
 - \em score (\c pyprophet score)\n
 - (optional) \em infer peptide/protein/gene (\c pyprophet infer peptide|protein|gene)\n
 - (optional) \em infer peptidoform (\c pyprophet infer peptidoform)\n
 - (optional) \em export reports/plots/Proteomics TSV/Proteomics matrix/Small-molecule TSV (\c pyprophet export ...)\n

It requires PyProphet standalone executables or a \c pyprophet on the PATH.

Multithreading: the global \c -threads parameter is passed through to PyProphet.

\par Merging behavior
If multiple \c -in OSW files are provided, a merge step is executed first via \c pyprophet merge osw.\n
Inputs are passed positionally as \c [INFILES]..., with the first input used as the \c --template OSW.\n
If a single input is provided and \c -out is omitted, scoring is performed in place (the input file is updated).

\par Scoring (semi-supervised)
Supports multiple levels (\c ms1, \c ms2, \c ms1ms2, \c transition, \c alignment) and classifiers (\c LDA, \c SVM, \c XGBoost, \c HistGradientBoosting).\n
Additional options include subsampling, cross-validation/semi-supervised iterations, initial/iteration FDRs, optional main score override, feature scaling, and an \c --autotune switch for supported classifiers.\n
Pre-trained weights can be applied via \c scoring:apply_weights to skip training.

\par Exports
\li Reports (single-file PDF) and score plots (PDF) do not require an explicit output path.\n
\li Proteomics TSV: enable with \c export:run_tsv; the output path is &lt;out&gt;.tsv.\n
\li Proteomics matrix: enable with \c export:run_matrix; the output path is &lt;out&gt;.matrix.tsv.\n
\li Small molecules TSV: enable with \c export:run_compound; the output path is &lt;out&gt;.compound.tsv.\n
    Format selection: \c export:compound:format = \c matrix or \c legacy_merged.\n
    Filtering: \c export:compound:max_rs_peakgroup_qvalue limits by run-specific peak group-level q-value.\n
\li Mutual exclusion: \c export:run_compound cannot be combined with \c export:run_tsv or \c export:run_matrix.

\par Logging and dry runs
The adapter captures and echoes PyProphet stdout/stderr so loguru messages are visible at default verbosity.\n
With \c -dry_run, the constructed PyProphet command(s) are printed but not executed.

\par Note on HistGradientBoosting threading
When using \c --classifier \c HistGradientBoosting for scoring, the \c OMP_NUM_THREADS environment variable controls OpenMP thread usage to avoid CPU oversubscription.\n
PyProphet may set a default if not specified, but for best control/performance set it explicitly before launching PyProphet.\n
Example: on a machine with 20 CPU threads and \c -threads \c 3 for semi-supervised learning, set \c OMP_NUM_THREADS=6 (floor(20/3)).

\par Note on Parquet support
PyProphet can read/write Parquet; this adapter currently supports only SQLite-based OSW files.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PyProphetAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PyProphetAdapter.html
*/



/// @cond TOPPCLASSES
using namespace OpenMS;
using std::vector;

class PyProphetAdapter final : public TOPPBase
{
public:
  PyProphetAdapter() :
      TOPPBase("PyProphetAdapter", "Run PyProphet (score/IPF/export) on OSW results.", true,
               {
                 {"Roest, H.L. et al.",
                  "OpenSWATH enables automated, targeted analysis of data-independent acquisition MS data",
                  "Nature Biotechnology volume 32, pages 219–223 (2014)",
                  "https://doi.org/10.1038/nbt.2841"},
                 {"Rosenberger, G. et al.",
                  "Inference and quantification of peptidoforms in large sample cohorts by SWATH-MS",
                  "Nature Biotechnology volume 35, pages 781–788 (2017)",
                  "https://doi.org/10.1038/nbt.3908"},
                 {"Meier, F. et al.",
                  "diaPASEF: parallel accumulation–serial fragmentation combined with data-independent acquisition",
                  "Nature Methods volume 17, pages 1229–1236 (2020)",
                  "https://doi.org/10.1038/s41592-020-00998-0"}
               })
  {}

protected:
  // --- Scoring --------------------------------------
  // Classifier enum for supported models in PyProphet
  enum class Classifier
  {
    LDA = 0,
    SVM,
    XGBoost,
    HistGradientBoosting
  };

  // canonical strings for CLI + INI validation (order must match enum)
  static const StringList& classifier_names_()
  {
    static StringList names = ListUtils::create<String>("LDA,SVM,XGBoost,HistGradientBoosting");
    return names;
  }

  // convert enum -> canonical string
  static String classifier_to_string_(Classifier c)
  {
    const auto& names = classifier_names_();
    switch (c)
    {
      case Classifier::LDA:                  return names[0];
      case Classifier::SVM:                  return names[1];
      case Classifier::XGBoost:              return names[2];
      case Classifier::HistGradientBoosting: return names[3];
    }
    // Fallback
    return names[1];
  }

  // parse user string (case-insensitive) -> enum
  static bool parse_classifier_(const String& s, Classifier& out)
  {
    const String ss = String(s).trim().toLower();
    if (ss == "lda") { out = Classifier::LDA; return true; }
    if (ss == "svm") { out = Classifier::SVM; return true; }
    if (ss == "xgboost" || ss == "xgb") { out = Classifier::XGBoost; return true; }
    if (ss == "histgradientboosting" || ss == "histgb" || ss == "hgb") { out = Classifier::HistGradientBoosting; return true; }
    return false;
  }

  // Levels for PyProphet scoring
  enum class Level
  {
    MS1 = 0,
    MS2,
    MS1MS2,
    TRANSITION,
    ALIGNMENT
  };

  static const StringList& level_names_()
  {
    // canonical, lowercase (what PyProphet expects)
    static StringList names = ListUtils::create<String>("ms1,ms2,ms1ms2,transition,alignment");
    return names;
  }

  static String level_to_string_(Level L)
  {
    const auto& n = level_names_();
    switch (L)
    {
      case Level::MS1:        return n[0];
      case Level::MS2:        return n[1];
      case Level::MS1MS2:     return n[2];
      case Level::TRANSITION: return n[3];
      case Level::ALIGNMENT:  return n[4];
    }
    // Fallback to ms2
    return n[1];
  }

  static bool parse_level_token_(const String& tok, Level& out)
  {
    String s = String(tok).trim().toLower();
    if (s == "ms1")        { out = Level::MS1;        return true; }
    if (s == "ms2")        { out = Level::MS2;        return true; }
    if (s == "ms1ms2")     { out = Level::MS1MS2;     return true; }
    if (s == "transition") { out = Level::TRANSITION; return true; }
    if (s == "alignment")  { out = Level::ALIGNMENT;  return true; }
    return false;
  }

  /**
   * Parse CSV 'level' string, validate combo, and return a canonical ordered list.
   *
   * Rules enforced:
   *  - 'ms1ms2' cannot be combined with 'ms2' in the same run.
   *  - 'transition' require at least one base level ('ms1'|'ms2'|'ms1ms2').
   *  - Duplicates are removed; final order is canonical: ms1, ms2, ms1ms2, transition, alignment.
   */
  static std::vector<Level> normalize_and_validate_levels_(const String& csv, String& error)
  {
    // parse tokens (CSV, whitespace-insensitive)
    StringList tokens; csv.split(',', tokens);

    bool seen_ms1=false, seen_ms2=false, seen_ms1ms2=false, seen_tr=false, seen_al=false;
    StringList bad;

    for (const auto& tok : tokens)
    {
      String t = tok;
      t.trim();
      if (t.trim().empty()) continue;
      Level L;
      if (!parse_level_token_(t, L))
      {
        bad.push_back(t);
        continue;
      }
      switch (L)
      {
        case Level::MS1:        seen_ms1 = true; break;
        case Level::MS2:        seen_ms2 = true; break;
        case Level::MS1MS2:     seen_ms1ms2 = true; break;
        case Level::TRANSITION: seen_tr = true; break;
        case Level::ALIGNMENT:  seen_al = true; break;
      }
    }

    if (!bad.empty())
    {
      error = "Unknown level(s): " + ListUtils::concatenate(bad, ",") +
              ". Allowed: " + ListUtils::concatenate(level_names_(), ",");
      return {};
    }

    // Default if nothing was provided (default is "ms2")
    if (!seen_ms1 && !seen_ms2 && !seen_ms1ms2 && !seen_tr && !seen_al)
    {
      seen_ms2 = true;
    }

    // Combo checks
    if (seen_ms1ms2 && seen_ms2)
    {
      error = "Invalid combination: 'ms1ms2' cannot be combined with 'ms2' in the same run.";
      return {};
    }

    const bool any_base = (seen_ms1 || seen_ms2 || seen_ms1ms2);
    if ((seen_tr) && !any_base)
    {
      error = "Invalid combination: 'transition' require at least one base level ('ms1', 'ms2', or 'ms1ms2').";
      return {};
    }

    // Canonical order for execution
    std::vector<Level> out;
    if (seen_ms1)    out.push_back(Level::MS1);
    if (seen_ms2)    out.push_back(Level::MS2);
    if (seen_ms1ms2) out.push_back(Level::MS1MS2);
    if (seen_tr)     out.push_back(Level::TRANSITION);
    if (seen_al)     out.push_back(Level::ALIGNMENT);
    return out;
  }

  // --- Inference Context --------------------------------------

  // Context to estimate peptide-level FDR control.
  enum class Context
  {
    RUN_SPECIFIC = 0,
    EXPERIMENT_WIDE,
    GLOBAL
  };

  static const StringList& context_names_()
  {
    // canonical values as expected by PyProphet CLI
    static StringList names = ListUtils::create<String>("run-specific,experiment-wide,global");
    return names;
  }

  static String context_to_string_(Context c)
  {
    const auto& n = context_names_();
    switch (c)
    {
      case Context::RUN_SPECIFIC:   return n[0];
      case Context::EXPERIMENT_WIDE:return n[1];
      case Context::GLOBAL:         return n[2];
    }
    return n[2];
  }

  static bool parse_context_token_(const String& tok, Context& out)
  {
    String s = String(tok).trim().toLower();
    if (s == "run-specific" || s == "runspecific")       { out = Context::RUN_SPECIFIC;    return true; }
    if (s == "experiment-wide" || s == "experimentwide") { out = Context::EXPERIMENT_WIDE; return true; }
    if (s == "global")                                    { out = Context::GLOBAL;          return true; }
    return false;
  }

  /**
   * Parse CSV 'context' string and return validated, order-preserving unique list.
   * Defaults to {"run-specific"} if empty.
   */
  static std::vector<Context> parse_and_validate_contexts_(const String& csv, String& error)
  {
    StringList toks; csv.split(',', toks);
    std::vector<Context> out;
    std::unordered_set<int> seen;

    StringList bad;
    for (const auto& tok : toks)
    {
      String t = tok; t.trim();
      if (t.empty()) continue;

      Context c;
      if (!parse_context_token_(t, c)) { bad.push_back(t); continue; }

      const int key = static_cast<int>(c);
      if (seen.insert(key).second) out.push_back(c);
    }

    if (!bad.empty())
    {
      error = "Unknown context(s): " + ListUtils::concatenate(bad, ",") +
              ". Allowed: " + ListUtils::concatenate(context_names_(), ",");
      return {};
    }

    if (out.empty()) out.push_back(Context::RUN_SPECIFIC);
    return out;
  }

  // helper for "true"/"false" string options
  static bool is_true_(const String& v)
  {
    String s = v; s.trim(); s.toLower();
    return s == "true";
  }

  // ---- helpers for argument handling and dry-run printing ----

  // Parse a shell-like argument string into tokens (supports quotes and backslash escapes)
  // Additionally supports escaping leading dashes to get past the OpenMS CLI parser:
  //   - Prefix a leading dash with '+'  :  "+--flag"  -> "--flag"
  //   - Or escape with backslash        :  "\-flag"   -> "-flag"
  // You can use multiple tokens inside quotes, e.g. "+--a 1 +--b 2"
  static vector<String> tokenize_args_(const String& extra)
  {
    vector<String> out;
    const std::string s = extra;
    std::string cur;
    bool in_single = false, in_double = false, escaping = false;

    for (char ch : s)
    {
      if (escaping) { cur += ch; escaping = false; continue; }
      if (ch == '\\') { escaping = true; continue; }
      if (!in_double && ch == '\'') { in_single = !in_single; continue; }
      if (!in_single && ch == '\"') { in_double = !in_double; continue; }
      if (!in_single && !in_double && std::isspace(static_cast<unsigned char>(ch)))
      {
        if (!cur.empty()) { out.emplace_back(cur); cur.clear(); }
        continue;
      }
      cur += ch;
    }
    if (!cur.empty()) out.emplace_back(cur);

    // Unescape leading dashes on each token:
    // "+--foo" -> "--foo", "+-foo" -> "-foo", "\-foo" -> "-foo"
    for (auto& tok : out)
    {
      std::string t = tok;
      if (t.rfind("+--", 0) == 0)      { t.erase(0, 1); }   // drop the '+'
      else if (t.rfind("+-", 0) == 0)  { t.erase(0, 1); }   // drop the '+'
      else if (t.rfind("\\-", 0) == 0) { t.erase(0, 1); }   // drop the backslash
      tok = t;
    }
    return out;
  }

  static String shell_quote_(const String& s)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    bool need = s.find_first_of(" \t\"") != String::npos;
    if (!need) return s;
    String r = s; r.substitute("\"", "\\\"");
    return "\"" + r + "\"";
#else
    bool need = s.find_first_of(" \t\"'\\$`") != String::npos;
    if (!need) return s;
    String r = s; r.substitute("'", "'\"'\"'");
    return "'" + r + "'";
#endif
  }

  static String join_cmd_(const String& exe, const vector<String>& args)
  {
    String cmd = shell_quote_(exe);
    for (const auto& a : args) { cmd += " "; cmd += shell_quote_(a); }
    return cmd;
  }

  void registerOptionsAndFlags_() override
  {
    // Inputs: support one or more OSW (if >1, we'll merge first)
    registerInputFileList_("in", "<files>", StringList(),
                           "Input OSW file(s). If multiple are provided, a merge step is executed first.", true);
    setValidFormats_("in", ListUtils::create<String>("OSW"));

    // Output OSW
    registerOutputFile_("out", "<file>", "",
                        "Output OSW file. If empty and only one input is provided, scoring is done in place.", false);
    setValidFormats_("out", ListUtils::create<String>("OSW"));

    // Scoring
    registerTOPPSubsection_("scoring", "Scoring module of PyProphet for semi-supervised learning");
    registerStringOption_("scoring:run_score", "<true|false>", "true","Run semi-supervised scoring.", false, false);
    setValidStrings_("scoring:run_score", ListUtils::create<String>("true,false"));
    registerStringOption_("scoring:level", "<level[,level,...]>", "ms1ms2",
                          "OSW data level(s) for scoring. Accepts a single level (e.g., 'ms2') or a comma-separated list (e.g., 'ms1,ms2,transition'). "
                          "Allowed values: ms1, ms2, ms1ms2, transition, alignment.", false);

    registerStringOption_("scoring:classifier", "<name>", "SVM",
                          "Classifier for semi-supervised scoring.", false);
    setValidStrings_("scoring:classifier", classifier_names_());

    registerDoubleOption_("scoring:subsample_ratio", "<0..1>", 1.0,
                          "Subsample ratio for training on large datasets; weights are applied to the full set.", false, true);
    setMinFloat_("scoring:subsample_ratio", 0.0);
    setMaxFloat_("scoring:subsample_ratio", 1.0);

    registerInputFile_("scoring:apply_weights", "<file>", "",
                       "Apply a pre-trained weights file (skip training). The same file will be applied to all requested levels.", false);
    registerFlag_("scoring:autotune", "Autotune hyperparameters for XGBoost/SVM/HistGradientBoosting.", false);

    registerIntOption_("scoring:xeval_num_iter", "<number>", 10, "Number of iterations for cross-validation of semi-supervised learning step.", false, true);
    setMinInt_("scoring:xeval_num_iter", 1);
    registerIntOption_("scoring:ss_num_iter", "<number>", 10, "Number of iterations for semi-supervised learning step.", false, true);
    setMinInt_("scoring:ss_num_iter", 1);
    registerDoubleOption_("scoring:ss_initial_fdr", "<0..1>", 0.15,
                          "Initial FDR cutoff for best scoring targets.", false, true);
    setMinFloat_("scoring:ss_initial_fdr", 0.0);
    setMaxFloat_("scoring:ss_initial_fdr", 1.0);
    registerDoubleOption_("scoring:ss_iteration_fdr", "<0..1>", 0.05,
                          "Iteration FDR cutoff for best scoring targets.", false, true);
    setMinFloat_("scoring:ss_iteration_fdr", 0.0);
    setMaxFloat_("scoring:ss_iteration_fdr", 1.0);
    registerStringOption_("scoring:ss_main_score", "<var>", "auto",
                          "Main score to start semi-supervised-learning. Default is set to auto, meaning each iteration of\n"
                          "learning a dynamic main score selection process will occur. If you want to have a set starting\n"
                          "main score for each learning iteration, you can set a specifc score, i.e. \"var_xcorr_shape\"",
                          false,
                          true);
    registerFlag_("scoring:ss_scale_features", "Scale / standardize features to unit variance before semi-supervised learning.", false);

    registerStringOption_("scoring:extra", "<args>", "",
                          "Advanced: extra args passed verbatim to 'pyprophet score' (e.g. \"+--parametric +--pi0_lambda 0.4 0 0\"), note the `+` is needed infront of the `--` to avoid it being interpreted as an OpenMS TOPP tool parameter."
                          "Examples: -scoring:extra \"+--parametric +--pi0_lambda 0.4 0 0\"",
                          false,
                          true);

    // Peptide /Protein/ Gene Inference
    registerTOPPSubsection_("infer", "Inference module of PyProphet for error-rate control for Peptide/Peptidoform(IPF)/Protein/Gene levels");
    registerTOPPSubsection_("infer:peptide", "Error-rate control for peptide-level");
    registerStringOption_("infer:run_peptide", "<true|false>", "true","Run peptide inference before IPF.", false, false);
    setValidStrings_("infer:run_peptide", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:peptide:extra", "<args>", "",
                          "Advanced: extra args passed verbatim to 'pyprophet infer peptide'.",
                          false, true);

    registerTOPPSubsection_("infer:protein", "Error-rate control for protein-level");
    registerStringOption_("infer:run_protein", "<true|false>", "true","Run protein inference before IPF.", false, false);
    setValidStrings_("infer:run_protein", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:protein:extra", "<args>", "",
                          "Advanced: extra args passed verbatim to 'pyprophet infer protein'.",
                          false, true);

    registerTOPPSubsection_("infer:gene", "Error-rate control for gene-level");
    registerStringOption_("infer:run_gene", "<true|false>", "false","Run gene inference before IPF.", false, true);
    setValidStrings_("infer:run_gene", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:gene:extra", "<args>", "",
                          "Advanced: extra args passed verbatim to 'pyprophet infer gene'.",
                          false, true);

    // Inference context
    registerStringOption_("infer:context", "<context[,context,...]>", "global",
                          "Context(s) for peptide/protein/gene inference. Accepts a single value or a comma-separated list. "
                          "Allowed values: run-specific, experiment-wide, global.", false, false);

    // IPF
    registerStringOption_("infer:run_ipf", "<true|false>", "false","Run peptidoform inference (IPF) after scoring.", false, false);
    setValidStrings_("infer:run_ipf", ListUtils::create<String>("true,false"));
    registerStringOption_("infer:ipf_ms1_scoring", "<true|false>", "false","Use MS1 precursor information in IPF scoring.", false, true);
    setValidStrings_("infer:ipf_ms1_scoring", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:ipf_ms2_scoring", "<true|false>", "false","Use MS2 precursor information in IPF scoring.", false, true);
    setValidStrings_("infer:ipf_ms2_scoring", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:ipf_h0", "<true|false>", "false","Include possibility that peak groups are not covered by peptidoform space.", false, true);
    setValidStrings_("infer:ipf_h0", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:ipf_grouped_fdr", "<true|false>", "false","Compute grouped FDR instead of pooled FDR.", false, true);
    setValidStrings_("infer:ipf_grouped_fdr", ListUtils::create<String>("true,false"));

    registerDoubleOption_("infer:ipf_max_precursor_pep", "<double>", 0.7, "Max PEP for scored precursors considered in IPF.", false, true);
    registerDoubleOption_("infer:ipf_max_peakgroup_pep", "<double>", 0.7, "Max PEP for scored peak-groups considered in IPF.", false, true);
    registerDoubleOption_("infer:ipf_max_precursor_peakgroup_pep", "<double>", 0.4, "Max integrated precursor-peakgroup PEP.", false, true);
    registerDoubleOption_("infer:ipf_max_transition_pep", "<double>", 0.6, "Max PEP for scored transitions considered in IPF.", false, true);

    // Optional exports
    registerTOPPSubsection_("export", "Export module of PyProphet for various exports");
    registerStringOption_("export:score_report", "<true|false>", "true",
                          "Generate a report that summarizes the results of your analysis, including scores and identifications, and other relevant information.", false, false);
    setValidStrings_("export:score_report", ListUtils::create<String>("true,false"));
    registerStringOption_("export:score_plots", "<true|false>", "false",
                          "Export score plots (PDF).", false, true);
    setValidStrings_("export:score_plots", ListUtils::create<String>("true,false"));

    registerTOPPSubsection_("export:tsv", "Export long-format TSV.");
    registerStringOption_("export:run_tsv", "<true|false>", "true",
                          "Export Proteomics TSV results (long format).", false, false);
    setValidStrings_("export:run_tsv", ListUtils::create<String>("true,false"));
    registerStringOption_("export:tsv:extra", "<args>", "",
                          "Advanced: extra args passed verbatim to 'pyprophet export tsv'.",
                          false, true);

    registerTOPPSubsection_("export:compound", "Export Small Molecules long-format TSV.");
    registerStringOption_("export:run_compound", "<true|false>", "false",
                          "Export  Small Molecules TSV results.", false, false);
    setValidStrings_("export:run_compound", ListUtils::create<String>("true,false"));
    registerStringOption_("export:compound:format", "<matrix|legacy_merged>", "legacy_merged",
                          "Export format, either matrix, legacy_merged (PyProphet) ", false, true);
    setValidStrings_("export:compound:format", ListUtils::create<String>("matrix,legacy_merged"));
    registerDoubleOption_("export:compound:max_rs_peakgroup_qvalue", "<double>", 0.05,
                          "Filter results to maximum run-specific peak group-level q-value.", false, true);


    registerTOPPSubsection_("export:matrix", "Export matrix TSV.");
    registerStringOption_("export:run_matrix", "<true|false>", "true",
                          "Export TSV quantification matrix.", false, false);
    setValidStrings_("export:run_matrix", ListUtils::create<String>("true,false"));
    registerStringOption_("export:matrix:extra", "<args>", "",
                          "Advanced: extra args passed verbatim to 'pyprophet export matrix'.",
                          false, true);


    // Executable
    registerInputFile_("pyprophet_executable", "<executable>",
#ifdef OPENMS_WINDOWSPLATFORM
                       "pyprophet.exe",
#else
                       "pyprophet",
#endif
                       "PyProphet executable (full/relative path or found via PATH).", false, false, {"is_executable"});

    // Logging / misc
    registerStringOption_("log_level", "<level>", "INFO",
                          "PyProphet log level (TRACE, DEBUG, INFO, SUCCESS, WARNING, ERROR, CRITICAL).", false);
    setValidStrings_("log_level", ListUtils::create<String>("TRACE,DEBUG,INFO,SUCCESS,WARNING,ERROR,CRITICAL"));

    registerFlag_("dry_run", "Print commands without executing.", false);
  }

  // ---- process runner (TOPPBase; capture & echo pyprophet logs) -------------
  ExitCodes run_pyprophet_(const String& exe, const std::vector<String>& args)
  {
    // Build Qt arg list for TOPPBase runner
    QStringList qargs;
    qargs.reserve(static_cast<int>(args.size()));
    for (const auto& a : args) qargs << a.toQString();

    // Capture process output
    String proc_out, proc_err;
    ExitCodes ec = runExternalProcess_(exe.toQString(), qargs, proc_out, proc_err);

    // Echo stdout (PyProphet/loguru often logs here)
    if (!proc_out.empty())
    {
      // print as INFO so users see it at default verbosity
      writeLogInfo_(proc_out);
    }

    // Echo stderr:
    // - if success, show as INFO (so users still see warnings/progress)
    // - if failure, show as ERROR
    if (!proc_err.empty())
    {
      if (ec == EXECUTION_OK) writeLogInfo_(proc_err);
      else                    writeLogError_(proc_err);
    }

    return ec;
  }

  ExitCodes main_(int, const char**) override
  {
    // Inputs
    StringList in_list = getStringList_("in");
    if (in_list.empty())
    {
      writeLogError_("No input (-in) given.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // Output
    String out = getStringOption_("out");
    if (out.empty() && in_list.size() == 1) out = in_list.front(); // in-place scoring by default, non-destructive
    if (out.empty())
    {
      writeLogError_("Multiple inputs provided but no -out given.");
      return ILLEGAL_PARAMETERS;
    }

    FileTypes::Type out_type = FileHandler::getTypeByFileName(out);
    if (out_type == FileTypes::UNKNOWN) out_type = FileTypes::OSW;
    if (out_type != FileTypes::OSW)
    {
      writeLogError_("Output must be OSW for PyProphetAdapter.");
      return PARSE_ERROR;
    }

    const String exe = getStringOption_("pyprophet_executable");

    const bool do_score        = is_true_(getStringOption_("scoring:run_score"));
    const String level_csv = getStringOption_("scoring:level");

    // parse classifier option into enum
    const String classifier_str = getStringOption_("scoring:classifier");
    Classifier classifier_enum{};
    if (!parse_classifier_(classifier_str, classifier_enum))
    {
      writeLogError_("Unknown classifier '" + classifier_str + "'. Allowed: " +
                     ListUtils::concatenate(classifier_names_(), ","));
      return ILLEGAL_PARAMETERS;
    }
    const String classifier_cli = classifier_to_string_(classifier_enum); // canonical CLI string

    const double subsample = getDoubleOption_("scoring:subsample_ratio");
    const String weights = getStringOption_("scoring:apply_weights");
    const bool   score_autotune         = getFlag_("scoring:autotune");
    const int    score_xeval_num_iter   = getIntOption_("scoring:xeval_num_iter");
    const int    score_ss_num_iter      = getIntOption_("scoring:ss_num_iter");
    const double score_ss_initial_fdr   = getDoubleOption_("scoring:ss_initial_fdr");
    const double score_ss_iteration_fdr = getDoubleOption_("scoring:ss_iteration_fdr");
    const String score_ss_main_score    = getStringOption_("scoring:ss_main_score");
    const bool   score_ss_scale_features= getFlag_("scoring:ss_scale_features");

    const bool do_ipf        = is_true_(getStringOption_("infer:run_ipf"));
    const bool do_infer_pep  = is_true_(getStringOption_("infer:run_peptide"));
    const bool do_infer_pro  = is_true_(getStringOption_("infer:run_protein"));
    const bool do_infer_gene = is_true_(getStringOption_("infer:run_gene"));

    const String log_level = getStringOption_("log_level");
    const bool dry_run = getFlag_("dry_run");
    const int threads = getIntOption_("threads");

    // Parse & validate levels (CSV -> enum vector)
    String error_msg;
    std::vector<Level> levels = normalize_and_validate_levels_(getStringOption_("scoring:level"), error_msg);
    if (!error_msg.empty())
    {
      writeLogError_(error_msg);
      return ILLEGAL_PARAMETERS;
    }

    // (optional) merge when multiple inputs  ----------------------------------
    String merged_or_input = in_list.front();
    if (in_list.size() > 1)
    {
      // Use the first input file as the required template OSW
      const String template_osw = in_list.front();

      String merged = File::absolutePath(File::getTempDirectory()) + "/" + File::getUniqueName() + "_merged.osw";
      vector<String> args;

      args.emplace_back("--log-level"); args.emplace_back(log_level);
      args.emplace_back("--no-log-colorize");
      args.emplace_back("merge"); args.emplace_back("osw");
      args.emplace_back("--out"); args.emplace_back(merged);
      args.emplace_back("--template"); args.emplace_back(template_osw);

      for (const auto& f : in_list)
      {
        args.emplace_back(f);
      }

      writeLogInfo_("====== Merging Input OSWs ===========================================");
      if (dry_run)
      {
        writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
      }
      else
      {
        ExitCodes ec = run_pyprophet_(exe, args);
        if (ec != EXECUTION_OK) return ec;
      }
      merged_or_input = merged;
    }

    //  score  -----------------------------------------------------------------
    //  run once per requested level, first pass reads merged/input; subsequent passes read 'out'
    if (do_score)
    {
      String current_in = merged_or_input;
      for (const Level lvl : levels)
      {
        const String lvl_cli = level_to_string_(lvl);

        vector<String> args;
        args.emplace_back("--log-level"); args.emplace_back(log_level);
        args.emplace_back("--no-log-colorize");
        args.emplace_back("score");
        args.emplace_back("--in"); args.emplace_back(current_in);
        args.emplace_back("--out"); args.emplace_back(out);
        args.emplace_back("--level"); args.emplace_back(lvl_cli);
        args.emplace_back("--classifier"); args.emplace_back(classifier_cli);
        args.emplace_back("--subsample_ratio"); args.emplace_back(String(subsample));
        args.emplace_back("--threads"); args.emplace_back(String(threads));

        if (!weights.empty()) { args.emplace_back("--apply_weights"); args.emplace_back(weights); }

        if (score_autotune) args.emplace_back("--autotune");
        args.emplace_back("--xeval_num_iter");   args.emplace_back(String(score_xeval_num_iter));
        args.emplace_back("--ss_num_iter");      args.emplace_back(String(score_ss_num_iter));
        args.emplace_back("--ss_initial_fdr");   args.emplace_back(String(score_ss_initial_fdr));
        args.emplace_back("--ss_iteration_fdr"); args.emplace_back(String(score_ss_iteration_fdr));

        if (!score_ss_main_score.empty())
        {
          args.emplace_back("--ss_main_score"); args.emplace_back(score_ss_main_score);
        }

        if (score_ss_scale_features)
        {
          args.emplace_back("--ss_scale_features");
        }

        {
          const vector<String> extra = tokenize_args_(getStringOption_("scoring:extra"));
          for (const auto& e : extra) args.emplace_back(e);
        }

        writeLogInfo_("====== Scoring level: " + lvl_cli + " ===========================================");
        if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
        else
        {
          ExitCodes ec = run_pyprophet_(exe, args);
          if (ec != EXECUTION_OK) return ec;
        }

        current_in = out;
      }
    }

    // Peptide / Protein / Gene inference --------------------------------------
    if (do_infer_pep || do_infer_pro || do_infer_gene)
    {
      String ctx_err;
      const auto contexts = parse_and_validate_contexts_(getStringOption_("infer:context"), ctx_err);
      if (!ctx_err.empty())
      {
        writeLogError_(ctx_err);
        return ILLEGAL_PARAMETERS;
      }

      auto maybe_run_infer = [&](const String& subject, const bool enabled, const String& ctx) -> ExitCodes
      {
        if (!enabled) return EXECUTION_OK;

        vector<String> args;
        args.emplace_back("--log-level"); args.emplace_back(log_level);
        args.emplace_back("--no-log-colorize");
        args.emplace_back("infer"); args.emplace_back(subject);
        args.emplace_back("--in");  args.emplace_back(out);
        args.emplace_back("--out"); args.emplace_back(out);
        args.emplace_back("--context"); args.emplace_back(ctx);

        // Append advanced passthrough args for this infer subcommand
        if (subject == "peptide")
        {
          const vector<String> extra = tokenize_args_(getStringOption_("infer:peptide:extra"));
          for (const auto& e : extra) args.emplace_back(e);
        }
        else if (subject == "protein")
        {
          const vector<String> extra = tokenize_args_(getStringOption_("infer:protein:extra"));
          for (const auto& e : extra) args.emplace_back(e);
        }
        else if (subject == "gene")
        {
          const vector<String> extra = tokenize_args_(getStringOption_("infer:gene:extra"));
          for (const auto& e : extra) args.emplace_back(e);
        }

        writeLogInfo_("====== " + subject + " Inference: " + ctx + " ===========================================");
        if (dry_run)
        {
          writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
          return EXECUTION_OK;
        }
        return run_pyprophet_(exe, args);
      };

      for (const auto c : contexts)
      {
        const String ctx_cli = context_to_string_(c);

        if (maybe_run_infer("peptide", do_infer_pep, ctx_cli) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
        if (maybe_run_infer("protein", do_infer_pro, ctx_cli) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
        if (maybe_run_infer("gene",    do_infer_gene, ctx_cli) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
      }
    }

    // Inference of Peptidoforms
    if (do_ipf)
    {
      const bool ipf_ms1    = is_true_(getStringOption_("infer:ipf_ms1_scoring"));
      const bool ipf_ms2    = is_true_(getStringOption_("infer:ipf_ms2_scoring"));
      const bool ipf_h0     = is_true_(getStringOption_("infer:ipf_h0"));
      const bool ipf_grouped= is_true_(getStringOption_("infer:ipf_grouped_fdr"));
      const double p1 = getDoubleOption_("infer:ipf_max_precursor_pep");
      const double p2 = getDoubleOption_("infer:ipf_max_peakgroup_pep");
      const double p3 = getDoubleOption_("infer:ipf_max_precursor_peakgroup_pep");
      const double p4 = getDoubleOption_("infer:ipf_max_transition_pep");

      vector<String> args;
      args.emplace_back("--log-level"); args.emplace_back(log_level);
      args.emplace_back("--no-log-colorize");
      args.emplace_back("infer"); args.emplace_back("peptidoform");
      args.emplace_back("--in");  args.emplace_back(out);
      args.emplace_back("--out"); args.emplace_back(out);
      args.emplace_back(ipf_ms1 ? "--ipf_ms1_scoring" : "--no-ipf_ms1_scoring");
      args.emplace_back(ipf_ms2 ? "--ipf_ms2_scoring" : "--no-ipf_ms2_scoring");
      args.emplace_back(ipf_h0  ? "--ipf_h0"          : "--no-ipf_h0");
      args.emplace_back(ipf_grouped ? "--ipf_grouped_fdr" : "--no-ipf_grouped_fdr");
      args.emplace_back("--ipf_max_precursor_pep");           args.emplace_back(String(p1));
      args.emplace_back("--ipf_max_peakgroup_pep");           args.emplace_back(String(p2));
      args.emplace_back("--ipf_max_precursor_peakgroup_pep"); args.emplace_back(String(p3));
      args.emplace_back("--ipf_max_transition_pep");          args.emplace_back(String(p4));

      writeLogInfo_("====== Peptidoform Inference ===========================================");
      if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
      else
      {
        ExitCodes ec = run_pyprophet_(exe, args);
        if (ec != EXECUTION_OK) return ec;
      }
    }

    // exporters --------------------------------------------------------------
    {
      const bool do_score_report = is_true_(getStringOption_("export:score_report"));
      const bool do_score_plots  = is_true_(getStringOption_("export:score_plots"));
      const bool do_tsv          = is_true_(getStringOption_("export:run_tsv"));
      const bool do_matrix       = is_true_(getStringOption_("export:run_matrix"));
      const bool do_compound     = is_true_(getStringOption_("export:run_compound"));

      // exclusivity: compound vs proteomics tsv/matrix
      if (do_compound && (do_tsv || do_matrix))
      {
        writeLogError_("Options 'export:run_compound' and ('export:run_tsv'/'export:run_matrix') are mutually exclusive. "
                       "Please enable either small-molecule export or proteomics export, not both.");
        return ILLEGAL_PARAMETERS;
      }

      auto run_export_flag = [&](const String& sub, bool enabled) -> ExitCodes
      {
        if (!enabled) return EXECUTION_OK;
        vector<String> args;
        args.emplace_back("--log-level"); args.emplace_back(log_level);
        args.emplace_back("--no-log-colorize");
        args.emplace_back("export"); args.emplace_back(sub);
        args.emplace_back("--in"); args.emplace_back(out);
        writeLogInfo_("====== Exporting " + sub + " ===========================================");
        if (dry_run) { writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args)); return EXECUTION_OK; }
        return run_pyprophet_(exe, args);
      };

      // score-report / score-plots (no --out)
      if (run_export_flag("score-report", do_score_report) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
      if (run_export_flag("score-plots",  do_score_plots)  != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;

      // Small molecules TSV (exclusive with proteomics TSV/matrix)
      if (do_compound)
      {
        const String fmt   = getStringOption_("export:compound:format"); // "matrix" or "legacy_merged"
        const double max_q = getDoubleOption_("export:compound:max_rs_peakgroup_qvalue");
        const String out_path = File::replaceExtension(File::appendSuffix(out, ".compound"), ".tsv");

        vector<String> args;
        args.emplace_back("--log-level"); args.emplace_back(log_level);
        args.emplace_back("--no-log-colorize");
        args.emplace_back("export");
        args.emplace_back("compound");
        args.emplace_back("--format"); args.emplace_back(fmt);
        args.emplace_back("--in");  args.emplace_back(out);
        args.emplace_back("--out"); args.emplace_back(out_path);
        args.emplace_back("--max_rs_peakgroup_qvalue"); args.emplace_back(max_q);

        writeLogInfo_("====== Exporting Small Molecules TSV (" + fmt + ") ===========================================");
        if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
        else { if (run_pyprophet_(exe, args) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR; }
      }

      // Proteomics TSV (derive "<out>.tsv")
      if (do_tsv)
      {
        const String tsv_path = File::replaceExtension(out, ".tsv");
        vector<String> args;
        args.emplace_back("--log-level"); args.emplace_back(log_level);
        args.emplace_back("--no-log-colorize");
        args.emplace_back("export"); args.emplace_back("tsv");
        args.emplace_back("--in");  args.emplace_back(out);
        args.emplace_back("--out"); args.emplace_back(tsv_path);

        // optional passthrough args
        if (!getStringOption_("export:tsv:extra").empty())
        {
          const vector<String> extra = tokenize_args_(getStringOption_("export:tsv:extra"));
          for (const auto& e : extra) args.emplace_back(e);
        }

        writeLogInfo_("====== Exporting Proteomics TSV ===========================================");
        if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
        else { if (run_pyprophet_(exe, args) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR; }
      }

      // Proteomics matrix (derive "<out>.matrix.tsv")
      if (do_matrix)
      {
        const String matrix_path = File::replaceExtension(File::appendSuffix(out, ".matrix"), ".tsv");
        vector<String> args;
        args.emplace_back("--log-level"); args.emplace_back(log_level);
        args.emplace_back("--no-log-colorize");
        args.emplace_back("export"); args.emplace_back("matrix");
        args.emplace_back("--in");  args.emplace_back(out);
        args.emplace_back("--out"); args.emplace_back(matrix_path);

        if (!getStringOption_("export:matrix:extra").empty())
        {
          const vector<String> extra = tokenize_args_(getStringOption_("export:matrix:extra"));
          for (const auto& e : extra) args.emplace_back(e);
        }

        writeLogInfo_("====== Exporting Matrix TSV ===========================================");
        if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + join_cmd_(exe, args));
        else { if (run_pyprophet_(exe, args) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR; }
      }
    }

    writeLogInfo_("PyProphetAdapter finished successfully.");
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  PyProphetAdapter tool;
  return tool.main(argc, argv);
}
/// @endcond
