// Copyright ...
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
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QStringList>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>


#include <array>
#include <unordered_set>

/**
@page TOPP_PyProphetAdapter PyProphetAdapter
@brief Adapter to run PyProphet scoring, optional, peptide, protein, and peptidoform inference (IPF), and exports on OSW files.

PyProphetAdapter orchestrates:\n
  1. (optional) @em merge of multiple OSW inputs (`pyprophet merge osw`)\n
  2. @em score (`pyprophet score`)\n
  3. (optional) @em infer peptidoform (`pyprophet infer peptidoform`)\n
  4. (optional) @em export reports/plots/TSV/matrix (`pyprophet export ...`)

It requires PyProphet standalone executables or a `pyprophet` on the PATH.

Multithreading: the global `-threads` parameter is passed through to PyProphet.

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

  // --- Inference Contenxt --------------------------------------

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


  void registerOptionsAndFlags_() override
  {
    // Inputs: support one or more OSW (if >1, we'll merge first)
    registerInputFileList_("in", "<files>", StringList(), "Input OSW file(s). If multiple are provided, a merge step is executed first.", false);
    setValidFormats_("in", ListUtils::create<String>("OSW"));

    // Output OSW
    registerOutputFile_("out", "<file>", "", "Output OSW file. If empty and only one input is provided, scoring is done in place.", false);
    setValidFormats_("out", ListUtils::create<String>("OSW"));

    // Scoring
    registerStringOption_("score:level", "<level[,level,...]>", "ms2",
                          "OSW data level(s) for scoring. Accepts a single level (e.g., 'ms2') or a comma-separated list (e.g., 'ms1,ms2,transition'). "
                          "Allowed values: ms1, ms2, ms1ms2, transition, alignment.", false);

    registerStringOption_("score:classifier", "<name>", "LDA",
                          "Classifier for semi-supervised scoring.", false);
    setValidStrings_("score:classifier", classifier_names_());

    registerDoubleOption_("score:subsample_ratio", "<0..1>", 1.0, "Subsample ratio for training on large datasets; weights are applied to the full set.", false, true);
    setMinFloat_("score:subsample_ratio", 0.0);

    registerInputFile_("score:apply_weights", "<file>", "", "Apply a pre-trained weights file (skip training). The same file will be applied to all requested levels.", true);

    // Peptide /Protein/ Gene Inference
    registerStringOption_("infer:peptide", "<true|false>", "true","Run peptide inference before IPF.", false, true);
    setValidStrings_("infer:peptide", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:protein", "<true|false>", "true","Run protein inference before IPF.", false, true);
    setValidStrings_("infer:protein", ListUtils::create<String>("true,false"));

    registerStringOption_("infer:gene", "<true|false>", "false","Run gene inference before IPF.", false, true);
    setValidStrings_("infer:gene", ListUtils::create<String>("true,false"));

    // Inference context (already added previously)
    registerStringOption_("infer:context", "<context[,context,...]>", "run-specific",
                          "Context(s) for peptide/protein/gene inference. Accepts a single value or a comma-separated list. "
                          "Allowed values: run-specific, experiment-wide, global.", false);

    // IPF (optional)
    registerStringOption_("ipf", "<true|false>", "false","Run peptidoform inference after scoring.", false, true);
    setValidStrings_("ipf", ListUtils::create<String>("true,false"));
    registerStringOption_("ipf:ms1_scoring", "<true|false>", "false","Use MS1 precursor information in IPF scoring.", false, true);
    setValidStrings_("ipf:ms1_scoring", ListUtils::create<String>("true,false"));

    registerStringOption_("ipf:ms2_scoring", "<true|false>", "false","Use MS2 precursor information in IPF scoring.", false, true);
    setValidStrings_("ipf:ms2_scoring", ListUtils::create<String>("true,false"));

    registerStringOption_("ipf:h0", "<true|false>", "false","Include possibility that peak groups are not covered by peptidoform space.", false, true);
    setValidStrings_("ipf:h0", ListUtils::create<String>("true,false"));

    registerStringOption_("ipf:grouped_fdr", "<true|false>", "false","Compute grouped FDR instead of pooled FDR.", false, true);
    setValidStrings_("ipf:grouped_fdr", ListUtils::create<String>("true,false"));

    registerDoubleOption_("ipf:max_precursor_pep", "<val>", 0.7, "Max PEP for scored precursors considered in IPF.", true);
    registerDoubleOption_("ipf:max_peakgroup_pep", "<val>", 0.7, "Max PEP for scored peak-groups considered in IPF.", true);
    registerDoubleOption_("ipf:max_precursor_peakgroup_pep", "<val>", 0.4, "Max integrated precursor-peakgroup PEP.", true);
    registerDoubleOption_("ipf:max_transition_pep", "<val>", 0.6, "Max PEP for scored transitions considered in IPF.", true);

    // Optional exports
    registerOutputFile_("export:score_report", "<pdf>", "", "Export a single-file PDF score report.", true);
    setValidFormats_("export:score_report", ListUtils::create<String>("pdf"));
    registerOutputFile_("export:score_plots", "<pdf>", "", "Export score plots (PDF).", true);
    setValidFormats_("export:score_plots", ListUtils::create<String>("pdf"));
    registerOutputFile_("export:tsv", "<tsv>", "", "Export TSV results.", true);
    setValidFormats_("export:tsv", ListUtils::create<String>("tsv"));
    registerOutputFile_("export:matrix", "<tsv>", "", "Export TSV quantification matrix.", true);
    setValidFormats_("export:matrix", ListUtils::create<String>("tsv"));

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

    registerFlag_("dry_run", "Print commands without executing.", true);
  }

  ExitCodes run_pyprophet_(const String& exe, const QStringList& args)
  {
    return runExternalProcess_(exe.toQString(), args);
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
    if (out.empty() && in_list.size() == 1) out = in_list.front(); // in-place scoring by default
    if (out.empty())
    {
      writeLogError_("Multiple inputs provided but no -out given.");
      return ILLEGAL_PARAMETERS;
    }

    // Enforce OSW -> OSW
    FileTypes::Type out_type = FileHandler::getTypeByFileName(out);
    if (out_type == FileTypes::UNKNOWN) out_type = FileTypes::OSW;
    if (out_type != FileTypes::OSW)
    {
      writeLogError_("Output must be OSW for PyProphetAdapter.");
      return PARSE_ERROR;
    }

    const String exe = getStringOption_("pyprophet_executable");
    const String level_csv = getStringOption_("score:level");

    // NEW: parse classifier option into enum
    const String classifier_str = getStringOption_("score:classifier");
    Classifier classifier_enum{};
    if (!parse_classifier_(classifier_str, classifier_enum))
    {
      writeLogError_("Unknown classifier '" + classifier_str + "'. Allowed: " +
                     ListUtils::concatenate(classifier_names_(), ","));
      return ILLEGAL_PARAMETERS;
    }
    const String classifier_cli = classifier_to_string_(classifier_enum); // canonical CLI string

    const double subsample = getDoubleOption_("score:subsample_ratio");
    const String weights = getStringOption_("score:apply_weights");
    const String log_level = getStringOption_("log_level");
    const bool do_ipf        = is_true_(getStringOption_("ipf"));
    const bool do_infer_pep  = is_true_(getStringOption_("infer:peptide"));
    const bool do_infer_pro  = is_true_(getStringOption_("infer:protein"));
    const bool do_infer_gene = is_true_(getStringOption_("infer:gene"));

    const bool dry_run = getFlag_("dry_run");
    const int threads = getIntOption_("threads"); // provided by TOPPBase common option

    // Parse & validate levels (CSV -> enum vector)
    String error_msg;
    std::vector<Level> levels = normalize_and_validate_levels_(getStringOption_("score:level"), error_msg);
    if (!error_msg.empty())
    {
      writeLogError_(error_msg);
      return ILLEGAL_PARAMETERS;
    }

    // 0) (optional) merge when multiple inputs
    String merged_or_input = in_list.front();
    if (in_list.size() > 1)
    {
      // Use the first input file as the required template OSW
      const String template_osw = in_list.front();

      String merged = File::absolutePath(File::getTempDirectory()) + "/" + File::getUniqueName() + "_merged.osw";
      QStringList args;

      args << "--log-level" << log_level.toQString() << "--no-log-colorize"
           << "merge" << "osw"
           << "--out" << merged.toQString()
           << "--template" << template_osw.toQString();

      for (const auto& f : in_list) args << "--in" << f.toQString();

      if (dry_run)
      {
        writeLogInfo_("DRY-RUN: " + exe + " " + String(args.join(' ').toStdString()));
      }
      else
      {
        ExitCodes ec = run_pyprophet_(exe, args);
        if (ec != EXECUTION_OK) return ec;
      }
      merged_or_input = merged;
    }

    //  score — run once per requested level, first pass reads merged/input; subsequent passes read 'out'
    {
      String current_in = merged_or_input;
      for (const Level lvl : levels)
      {
        const String lvl_cli = level_to_string_(lvl);

        QStringList args;
        args << "--log-level" << log_level.toQString() << "--no-log-colorize" << "score"
             << "--in" << current_in.toQString()
             << "--out" << out.toQString()
             << "--level" << lvl_cli.toQString()
             << "--classifier" << classifier_cli.toQString()  // from your classifier enum logic
             << "--subsample_ratio" << QString::number(subsample)
             << "--threads" << QString::number(threads);

        if (!weights.empty()) args << "--apply_weights" << weights.toQString();

        writeLogInfo_("Scoring level: " + lvl_cli);
        if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + String(args.join(' ').toStdString()));
        else
        {
          ExitCodes ec = run_pyprophet_(exe, args);
          if (ec != EXECUTION_OK) return ec;
        }

        current_in = out;
      }
    }

    // Peptide / Protein / Gene inference ---------------------
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

        QStringList args;
        args << "--log-level" << log_level.toQString() << "--no-log-colorize" << "infer" << subject.toQString()
             << "--in"  << out.toQString()
             << "--out" << out.toQString()
             << "--context" << ctx.toQString();

        if (dry_run)
        {
          writeLogInfo_("DRY-RUN: " + exe + " " + String(args.join(' ').toStdString()));
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
      const bool ipf_ms1    = is_true_(getStringOption_("ipf:ms1_scoring"));
      const bool ipf_ms2    = is_true_(getStringOption_("ipf:ms2_scoring"));
      const bool ipf_h0     = is_true_(getStringOption_("ipf:h0"));
      const bool ipf_grouped= is_true_(getStringOption_("ipf:grouped_fdr"));
      const double p1 = getDoubleOption_("ipf_max_precursor_pep");
      const double p2 = getDoubleOption_("ipf_max_peakgroup_pep");
      const double p3 = getDoubleOption_("ipf_max_precursor_peakgroup_pep");
      const double p4 = getDoubleOption_("ipf_max_transition_pep");

      QStringList args;
      args << "--log-level" << log_level.toQString() << "--no-log-colorize"
           << "infer" << "peptidoform"
           << "--in" << out.toQString()
           << "--out" << out.toQString()
           << (ipf_ms1 ? "--ipf_ms1_scoring" : "--no-ipf_ms1_scoring")
           << (ipf_ms2 ? "--ipf_ms2_scoring" : "--no-ipf_ms2_scoring")
           << (ipf_h0  ? "--ipf_h0"          : "--no-ipf_h0")
           << (ipf_grouped ? "--ipf_grouped_fdr" : "--no-ipf_grouped_fdr")
           << "--ipf_max_precursor_pep" << QString::number(p1)
           << "--ipf_max_peakgroup_pep" << QString::number(p2)
           << "--ipf_max_precursor_peakgroup_pep" << QString::number(p3)
           << "--ipf_max_transition_pep" << QString::number(p4);

      if (dry_run) writeLogInfo_("DRY-RUN: " + exe + " " + String(args.join(' ').toStdString()));
      else
      {
        ExitCodes ec = run_pyprophet_(exe, args);
        if (ec != EXECUTION_OK) return ec;
      }
    }

    // exporters
    auto run_export = [&](const String& sub, const String& dest) -> ExitCodes
    {
      if (dest.empty()) return EXECUTION_OK;
      QStringList args;
      args << "--log-level" << log_level.toQString() << "--no-log-colorize"
           << "export" << sub.toQString()
           << "--in" << out.toQString()
           << "--out" << dest.toQString();
      if (dry_run) { writeLogInfo_("DRY-RUN: " + exe + " " + String(args.join(' ').toStdString())); return EXECUTION_OK; }
      return run_pyprophet_(exe, args);
    };

    if (run_export("score-report", getStringOption_("export:score_report")) != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
    if (run_export("score-plots",  getStringOption_("export:score_plots"))  != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
    if (run_export("tsv",          getStringOption_("export:tsv"))          != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;
    if (run_export("matrix",       getStringOption_("export:matrix"))       != EXECUTION_OK) return EXTERNAL_PROGRAM_ERROR;

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
