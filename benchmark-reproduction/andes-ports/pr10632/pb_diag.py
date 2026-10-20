"""Per-species bias and spread from the proteobench intermediate (min 3 runs): median signed epsilon
(observed - expected log2 A/B; ratio compression shows as a bias towards 0) and median |epsilon_precision|
(deviation from the species' empirical median, i.e. spread without the bias)."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import pb_score as s
from proteobench.modules.quant.benchmarking import run_benchmarking
print("file\tspecies\tn\tmedian_log2_A_vs_B\texpected\tmedian_eps(signed)\tmedian_abs_eps_precision")
for p in sys.argv[1:]:
    inter, _, _ = run_benchmarking(p, "quantms", s.USER_INPUT, None, s.MODULE_SETTINGS_DIRS[s.MODULE], s.MODULE,
                                   "precursor ion", default_cutoff_min_feature=3)
    d = inter[inter["nr_observed"] >= 3]
    for sp, g in d.groupby("species"):
        print(f"{p}\t{sp}\t{len(g)}\t{g['log2_A_vs_B'].median():.4f}\t{g['log2_expectedRatio'].iloc[0]:.4f}\t"
              f"{g['epsilon'].median():.4f}\t{g['epsilon_precision_median'].abs().median():.4f}", flush=True)
