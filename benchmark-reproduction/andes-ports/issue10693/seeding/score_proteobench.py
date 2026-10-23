"""Score an MSstats CSV (ProteomicsLFQ / MSstatsConverter output) with the official proteobench 0.18.6
DDA ion QExactive module ("quantms" format), min 3 runs, as in #10585.
Prints: file, precursors, median |eps|, E. coli median |eps|, median CV."""
import sys, warnings
warnings.filterwarnings("ignore")
from proteobench.modules.quant.benchmarking import run_benchmarking
from proteobench.modules.constants import MODULE_SETTINGS_DIRS
from proteobench.datapoint.quant_datapoint import QuantDatapointHYE

MODULE = "quant_lfq_DDA_ion_QExactive"
# descriptive only, does not affect the metrics
UI = {"software_name": "quantms", "software_version": "ProteomicsLFQ", "search_engine": "Comet",
      "search_engine_version": "2025.01", "ident_fdr_psm": 0.01, "ident_fdr_peptide": None, "ident_fdr_protein": 0.01,
      "enable_match_between_runs": True, "precursor_mass_tolerance": "10 ppm", "fragment_mass_tolerance": "0.02 Da",
      "enzyme": "Trypsin", "allowed_miscleavages": 2, "min_peptide_length": 7, "max_peptide_length": 40,
      "comments_for_plotting": ""}

def score(csv, min_runs=3):
    inter, _, _ = run_benchmarking(csv, "quantms", UI, None, MODULE_SETTINGS_DIRS[MODULE], MODULE, "precursor ion",
                                   default_cutoff_min_feature=min_runs)
    dp = QuantDatapointHYE.generate_datapoint(inter, "quantms", UI, default_cutoff_min_feature=min_runs)
    r = dp.results[min_runs]
    return (r.get("nr_prec", r.get("nr_feature")), r["median_abs_epsilon_global"], r["median_abs_epsilon_ECOLI"], r["CV_median"])

if __name__ == "__main__":
    for csv in sys.argv[1:]:
        n, eps, ec, cv = score(csv)
        print(f"{csv}\t{n}\t{eps:.4f}\t{ec:.4f}\t{cv:.4f}", flush=True)
