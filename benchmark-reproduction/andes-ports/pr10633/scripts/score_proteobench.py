"""Score an MSstats CSV (ProteomicsLFQ / MSstatsConverter) with the official
proteobench package (0.18.6), module quant_lfq_DDA_ion_QExactive, input format
'quantms'. Prints precursors, median |eps| (global, E. coli) and median CV at
min_nr_observed (default 3), as one TSV line."""
import argparse, collections, json, sys

import proteobench
from proteobench.modules.constants import MODULE_SETTINGS_DIRS
from proteobench.modules.quant.benchmarking import run_benchmarking

ap = argparse.ArgumentParser()
ap.add_argument("msstats_csv")
ap.add_argument("--label", default="")
ap.add_argument("--min-obs", type=int, default=3)
ap.add_argument("--header", action="store_true")
a = ap.parse_args()

module_id = "quant_lfq_DDA_ion_QExactive"
user_input = collections.defaultdict(str, {
    "software_version": "dev", "search_engine": "ProSE", "search_engine_version": "dev",
    "ident_fdr_psm": 0.01, "ident_fdr_peptide": 0, "ident_fdr_protein": 0.01,
    "enable_match_between_runs": True, "precursor_mass_tolerance": "20 ppm",
    "fragment_mass_tolerance": "20 ppm", "enzyme": "Trypsin/P", "allowed_miscleavages": 2,
    "min_peptide_length": 7, "max_peptide_length": 40,
})
intermediate, _, _ = run_benchmarking(
    input_file=a.msstats_csv, input_format="quantms", user_input=user_input, all_datapoints=None,
    parse_settings_dir=MODULE_SETTINGS_DIRS[module_id], module_id=module_id,
    precursor_column_name="precursor ion", default_cutoff_min_feature=a.min_obs)
from proteobench.datapoint.quant_datapoint import QuantDatapointHYE
m = QuantDatapointHYE.get_metrics(intermediate, a.min_obs)[a.min_obs]
cols = ["label", "proteobench", "min_obs", "precursors", "median_abs_eps", "median_abs_eps_ECOLI",
        "median_abs_eps_HUMAN", "median_abs_eps_YEAST", "CV_median"]
vals = [a.label, proteobench.__version__, a.min_obs, m["nr_feature"], m["median_abs_epsilon_global"],
        m.get("median_abs_epsilon_ECOLI"), m.get("median_abs_epsilon_HUMAN"), m.get("median_abs_epsilon_YEAST"),
        m["CV_median"]]
if a.header:
    print("\t".join(cols))
print("\t".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in vals))
