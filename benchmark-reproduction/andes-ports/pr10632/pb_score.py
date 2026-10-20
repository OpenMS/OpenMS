"""Score an MSstats CSV (ProteomicsLFQ -out_msstats / MSstatsConverter) with the official proteobench package
(0.18.6), DDA quant ion module QExactive, input format "quantms", min 3 runs. Prints one TSV row per file."""
import sys, warnings
warnings.filterwarnings("ignore")
from proteobench.modules.quant.benchmarking import run_benchmarking
from proteobench.modules.constants import MODULE_SETTINGS_DIRS
from proteobench.datapoint.quant_datapoint import QuantDatapointHYE

MODULE = "quant_lfq_DDA_ion_QExactive"
USER_INPUT = {"software_name": "quantms", "software_version": "OpenMS", "search_engine": "Comet",
              "search_engine_version": "2025.01", "ident_fdr_psm": 0.01, "ident_fdr_peptide": 0.0,
              "ident_fdr_protein": 0.01, "enable_match_between_runs": False,
              "precursor_mass_tolerance": "10 ppm", "fragment_mass_tolerance": "0.02 Da", "enzyme": "Trypsin",
              "allowed_miscleavages": 2, "min_peptide_length": 7, "max_peptide_length": 40}

def score(path, min_runs=3):
    inter, _, _ = run_benchmarking(path, "quantms", USER_INPUT, None, MODULE_SETTINGS_DIRS[MODULE], MODULE,
                                   "precursor ion", default_cutoff_min_feature=min_runs)
    dp = QuantDatapointHYE.generate_datapoint(inter, "quantms", USER_INPUT, default_cutoff_min_feature=min_runs)
    r = dp["results"][min_runs]
    n = r.get("nr_prec", r.get("nr_feature"))
    return [n, r["median_abs_epsilon_global"], r["median_abs_epsilon_ECOLI"], r["median_abs_epsilon_YEAST"],
            r["median_abs_epsilon_HUMAN"], r["median_abs_epsilon_eq_species"], r["CV_median"]]

if __name__ == "__main__":
    print("file\tprecursors\tmedian_abs_eps\tecoli\tyeast\thuman\teq_species\tCV_median")
    for p in sys.argv[1:]:
        print(p + "\t" + "\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in score(p)), flush=True)
