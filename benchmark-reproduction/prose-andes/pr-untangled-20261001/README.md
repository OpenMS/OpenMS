# Untangled ProSE PRs vs. develop — 20-file benchmark (2026-10-01)

Supplement to OpenMS issue [#10364](https://github.com/OpenMS/OpenMS/issues/10364). The stacked ProSE PRs
#10368, #10366 and #10365 were ported onto `develop` as standalone PRs **#10377**, **#10378** and **#10379**.
Each was benchmarked against `develop` on the frozen 20-file suite, together with **#10335 merged with current
develop**. Every arm is compared with the same `develop` base, never stacked on another PR.

## What was compared

| Build | Source | Notes |
| --- | --- | --- |
| `base` | `develop` `f5ea2d04` | Includes #10357 (final-window fix), #10350 (index rebuild) and #10358 (initial-M clipping, on by default). No #10335 scoring or deduplication. |
| `pr10368` | #10377 head `a44c46c4` | Standalone #10368 (positive-intensity cutoff). |
| `pr10366` | #10378 `e0d77a6a` (algorithms; head `241f9a70` only fixes a test for clang) | Standalone #10366 (self-trained ion priors). |
| `pr10365_scorer` | #10379 commit 1 `825c33bb` | Opt-in mass-accuracy scorer, ported from #10335 `48dc1c2`. |
| `pr10365` | #10379 head `16ff0344` | + #10365: fitted kernel and `calibration:enabled=auto` (new default). |
| `pr10335` | local merge of `codex/native-cid` `c1b2f202` with develop `f5ea2d04` | `pr10335-on-develop.patch`; merge commit `dffc0f72` is not pushed. Conflicts with #10357/#10350/#10358 resolved by keeping both features. |

All builds passed `HyperScore_test`, `FragmentIndex_test` and `ProSEAlgorithm_test` (+ `FragmentIonLikelihoodModel_test`
for #10378) in the harness before searching. `sources.json` and `build_records.json` give exact commits, trees and
binary hashes.

Arms (each changes only the listed settings relative to the frozen `native_dedup:<dataset>:default` parameters of
the reference, which set `calibration:enabled=false`):

| Arm | Build | Change | Files |
| --- | --- | --- | --- |
| `base` | base | – | 20 |
| `pr68` | pr10368 | – | 20 |
| `pr66_priors` | pr10366 | `annotate:self_trained_ion_priors=true` | 20 |
| `pr65_default` | pr10365 | `calibration:enabled=auto` (the PR's new default; scorer stays HyperScore) | 20 |
| `pr65s_mass` | pr10365_scorer | `scoring:method=mass_accuracy`, calibration off (fixed 7 ppm kernel) | 14 (ppm) |
| `pr65s_mass_cal` | pr10365_scorer | `mass_accuracy`, `calibration:enabled=true` (fixed kernel) | 14 |
| `pr65_mass_cal` | pr10365 | `mass_accuracy`, `calibration:enabled=true` (fitted kernel) | 14 |
| `pr10335` | pr10335 | – (the frozen parameters are #10335's own default configuration) | 20 |

`peptide:deduplicate=true` is part of the frozen parameters. Builds without #10335 do not define that parameter;
the driver skipped it and logged the skip (`ignored_parameters` in each `summary.json`). Calibration is off in the frozen
parameters, so `pr65_default` sets the new `auto` default explicitly.

## Protocol

Same as the 2026-09-30 ablation supplement:
- **Inputs:** 20 selected mzML files (8,000 DDA MS2 each) and the three target/decoy FASTAs. Every file and FASTA
  hash was checked before each search.
- **PIN export:** identical to the 2026-09-30 ablation. RT features are removed and native features kept.
- **Rescoring:** Percolator 3.09.0 with `-Y -U --seed 1/42/137`, run independently per file. Target PSMs are counted
  at q ≤ 0.01.
- **Aggregation:** mean over seeds within a file, then mean over files within an instrument group.
- **Native counts:** rank-one target–decoy competition, `(D+1)/T ≤ 0.01`.

**Harness:** algorithms are compiled from each tree and linked against pyOpenMS nightly `3.7.0.dev20261001`, built
from `e25ac61`. The only ProSE-relevant library change between `e25ac61` and `develop` is FragmentIndex, which is
compiled from source. Toolchain: GCC 13.3, Python 3.11, 4 OpenMP threads per search, two searches concurrently. Times
are not a speed benchmark.

### Main table: develop vs. each PR (defaults of each PR's opt-in feature where applicable)

| Instrument / acquisition | n | develop | #10377 (cutoff) | #10378 priors | #10379 default (auto calibration) | #10335 merged with develop |
| --- | --- | --- | --- | --- | --- | --- |
| Velos CID | 3 | 2600.3 | 2600.3 (+0.00%) | 2590.8 (-0.37%) | 2600.3 (+0.00%) | 2904.0 (+11.68%) |
| HF-X HCD | 3 | 4380.3 | 4380.3 (+0.00%) | 4393.6 (+0.30%) | 4234.2 (-3.34%) | 4390.4 (+0.23%) |
| Astral HCD | 3 | 2189.6 | 2189.6 (+0.00%) | 2290.3 (+4.60%) | 2301.1 (+5.09%) | 2186.3 (-0.15%) |
| Lumos HCD LFQ | 3 | 5496.1 | 5496.1 (+0.00%) | 5490.0 (-0.11%) | 5267.6 (-4.16%) | 5500.2 (+0.07%) |
| Lumos CID TMT | 3 | 2327.3 | 2327.3 (+0.00%) | 2348.8 (+0.92%) | 2327.3 (+0.00%) | 2398.1 (+3.04%) |
| Exploris 480 TMTpro | 3 | 2244.9 | 2244.9 (+0.00%) | 2284.4 (+1.76%) | 2258.6 (+0.61%) | 2251.4 (+0.29%) |
| timsTOF HT | 2 | 1015.0 | 1015.0 (+0.00%) | 1038.8 (+2.35%) | 1023.0 (+0.79%) | 1002.0 (-1.28%) |

### #10379 opt-in mass-accuracy scorer (high-resolution files only)

| Instrument / acquisition | n | develop | scorer, fixed 7 ppm, no calibration | scorer + calibration, fixed kernel | + fitted kernel (#10365) vs. fixed kernel | native TDC delta of the fitted kernel |
| --- | --- | --- | --- | --- | --- | --- |
| HF-X HCD | 3 | 4380.3 | 4395.2 (+0.34%) | 4248.1 (-3.02%) | 4247.0 (-0.03%) | -564.7 |
| Astral HCD | 3 | 2189.6 | 2279.9 (+4.13%) | 2284.6 (+4.34%) | 2257.4 (-1.19%) | -577.3 |
| Lumos HCD LFQ | 3 | 5496.1 | 5488.2 (-0.14%) | 5313.2 (-3.33%) | 5309.6 (-0.07%) | -161.3 |
| Exploris 480 TMTpro | 3 | 2244.9 | 2243.3 (-0.07%) | 2255.7 (+0.48%) | 2272.2 (+0.73%) | -135.0 |
| timsTOF HT | 2 | 1015.0 | 1021.3 (+0.62%) | 1026.5 (+1.13%) | 1027.7 (+0.11%) | +41.5 |

### Native rank-one TDC delta vs. develop (mean per file)

| Instrument / acquisition | #10377 | #10378 | #10379 default | #10335 |
| --- | --- | --- | --- | --- |
| Velos CID | +0.0 | +0.0 | +0.0 | +509.7 |
| HF-X HCD | +0.0 | +0.0 | -3.7 | +0.0 |
| Astral HCD | +0.0 | +0.0 | +331.0 | +0.0 |
| Lumos HCD LFQ | +0.0 | +0.0 | -142.3 | +0.0 |
| Lumos CID TMT | +0.0 | +0.0 | +0.0 | +245.3 |
| Exploris 480 TMTpro | +0.0 | +0.0 | +174.3 | +0.0 |
| timsTOF HT | +0.0 | +0.0 | -4.5 | +0.0 |

### Per file: Percolator seeds 1 / 42 / 137 (mean)

| Dataset | develop | #10377 | #10378 priors | #10379 auto cal. | scorer | scorer+cal. | fitted kernel | #10335 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `velos_125_R1` | 2480/2505/2602 (2529.0) | 2480/2505/2602 (2529.0) | 2442/2513/2616 (2523.7) | 2480/2505/2602 (2529.0) | – | – | – | 2939/2874/3012 (2941.7) |
| `velos_5000_R2` | 2718/2750/2708 (2725.3) | 2718/2750/2708 (2725.3) | 2677/2706/2694 (2692.3) | 2718/2750/2708 (2725.3) | – | – | – | 2974/2981/2950 (2968.3) |
| `velos_25000_R3` | 2529/2563/2548 (2546.7) | 2529/2563/2548 (2546.7) | 2560/2553/2556 (2556.3) | 2529/2563/2548 (2546.7) | – | – | – | 2786/2853/2767 (2802.0) |
| `hfx_A2` | 4425/4524/4483 (4477.3) | 4425/4524/4483 (4477.3) | 4458/4516/4499 (4491.0) | 4249/4252/4307 (4269.3) | 4459/4516/4495 (4490.0) | 4225/4268/4260 (4251.0) | 4228/4253/4275 (4252.0) | 4434/4502/4520 (4485.3) |
| `hfx_B1` | 4267/4320/4265 (4284.0) | 4267/4320/4265 (4284.0) | 4318/4376/4286 (4326.7) | 4217/4198/4226 (4213.7) | 4298/4341/4325 (4321.3) | 4219/4226/4192 (4212.3) | 4198/4180/4144 (4174.0) | 4326/4344/4260 (4310.0) |
| `hfx_B3` | 4373/4337/4429 (4379.7) | 4373/4337/4429 (4379.7) | 4380/4337/4372 (4363.0) | 4228/4224/4207 (4219.7) | 4379/4391/4353 (4374.3) | 4272/4273/4298 (4281.0) | 4335/4287/4323 (4315.0) | 4364/4332/4432 (4376.0) |
| `astral_A2` | 1910/1983/1973 (1955.3) | 1910/1983/1973 (1955.3) | 2101/2077/2079 (2085.7) | 2070/2101/2051 (2074.0) | 1951/2055/2042 (2016.0) | 2007/2002/2103 (2037.3) | 2110/1925/1943 (1992.7) | 2019/1963/1944 (1975.3) |
| `astral_B1` | 2412/2321/2386 (2373.0) | 2412/2321/2386 (2373.0) | 2452/2401/2418 (2423.7) | 2419/2414/2422 (2418.3) | 2358/2439/2338 (2378.3) | 2345/2468/2378 (2397.0) | 2336/2371/2389 (2365.3) | 2371/2280/2326 (2325.7) |
| `astral_B3` | 2247/2246/2228 (2240.3) | 2247/2246/2228 (2240.3) | 2385/2270/2430 (2361.7) | 2362/2438/2433 (2411.0) | 2466/2419/2451 (2445.3) | 2403/2390/2465 (2419.3) | 2440/2387/2416 (2414.3) | 2256/2204/2314 (2258.0) |
| `lumos_lfq_5192` | 5505/5484/5511 (5500.0) | 5505/5484/5511 (5500.0) | 5486/5476/5512 (5491.3) | 5091/5097/5091 (5093.0) | 5484/5481/5482 (5482.3) | 5211/5198/5189 (5199.3) | 5190/5202/5165 (5185.7) | 5501/5511/5500 (5504.0) |
| `lumos_lfq_5194` | 5481/5464/5465 (5470.0) | 5481/5464/5465 (5470.0) | 5456/5459/5464 (5459.7) | 5293/5293/5296 (5294.0) | 5465/5467/5446 (5459.3) | 5338/5326/5323 (5329.0) | 5317/5329/5328 (5324.7) | 5478/5482/5484 (5481.3) |
| `lumos_lfq_5199` | 5511/5515/5529 (5518.3) | 5511/5515/5529 (5518.3) | 5518/5504/5535 (5519.0) | 5414/5411/5422 (5415.7) | 5522/5511/5536 (5523.0) | 5395/5422/5417 (5411.3) | 5421/5418/5416 (5418.3) | 5513/5506/5527 (5515.3) |
| `lumos_tmt_5058` | 2215/2292/2248 (2251.7) | 2215/2292/2248 (2251.7) | 2251/2302/2288 (2280.3) | 2215/2292/2248 (2251.7) | – | – | – | 2290/2304/2261 (2285.0) |
| `lumos_tmt_5059` | 2289/2349/2376 (2338.0) | 2289/2349/2376 (2338.0) | 2309/2436/2349 (2364.7) | 2289/2349/2376 (2338.0) | – | – | – | 2375/2376/2411 (2387.3) |
| `lumos_tmt_5066` | 2386/2407/2384 (2392.3) | 2386/2407/2384 (2392.3) | 2408/2391/2405 (2401.3) | 2386/2407/2384 (2392.3) | – | – | – | 2528/2518/2520 (2522.0) |
| `eclipse_tmtpro_10855` | 2405/2407/2443 (2418.3) | 2405/2407/2443 (2418.3) | 2410/2453/2433 (2432.0) | 2384/2424/2415 (2407.7) | 2426/2385/2380 (2397.0) | 2404/2381/2440 (2408.3) | 2402/2424/2437 (2421.0) | 2426/2370/2403 (2399.7) |
| `eclipse_tmtpro_10858` | 2169/2124/2104 (2132.3) | 2169/2124/2104 (2132.3) | 2179/2155/2105 (2146.3) | 2186/2222/2193 (2200.3) | 2152/2137/2085 (2124.7) | 2170/2178/2218 (2188.7) | 2230/2219/2228 (2225.7) | 2133/2132/2154 (2139.7) |
| `eclipse_tmtpro_10863` | 2180/2152/2220 (2184.0) | 2180/2152/2220 (2184.0) | 2262/2301/2262 (2275.0) | 2181/2173/2149 (2167.7) | 2201/2208/2216 (2208.3) | 2174/2140/2196 (2170.0) | 2151/2147/2212 (2170.0) | 2233/2139/2273 (2215.0) |
| `tims_plasma_30min` | 1015/1050/1026 (1030.3) | 1015/1050/1026 (1030.3) | 1067/989/1058 (1038.0) | 1030/1028/1039 (1032.3) | 1054/1019/1026 (1033.0) | 985/1022/1017 (1008.0) | 994/1003/1028 (1008.3) | 1035/1019/1038 (1030.7) |
| `tims_plasma_50min` | 1028/991/980 (999.7) | 1028/991/980 (999.7) | 1055/1018/1046 (1039.7) | 1018/1015/1008 (1013.7) | 1040/971/1018 (1009.7) | 1024/1064/1047 (1045.0) | 1036/1058/1047 (1047.0) | 963/1018/939 (973.3) |


### Isolation checks (all passed)

- **#10377:** normalized PIN byte-identical to `develop` on all 20 files.
- **#10378:** native candidate/score TSV identical to `develop` on all 20 files. The PIN is byte-identical after
  removing only the three `ion_prior_*` columns, so the gains come from the new Percolator features.
- **Kernel pair:** `pr65s_mass_cal` and `pr65_mass_cal` have identical calibrated precursor and fragment tolerances on
  all 14 files, so their difference isolates the fitted kernel.
- `checks.tsv` lists every check.

### Intensity-scale control for #10377

The same derived inputs as on 2026-09-30 (byte-identical) were used: each spectrum is scaled by a power of two so its
base peak lies in [1,2). This preserves m/z, RT, precursors and scan IDs.

| File | Arm | Unscaled seeds | Scaled seeds | Identical PIN |
| --- | --- | --- | --- | --- |
| velos_125_R1 | develop | [2480, 2505, 2602] | [2243, 2049, 2192] | False |
| velos_125_R1 | #10377 | [2480, 2505, 2602] | [2480, 2505, 2602] | True |
| eclipse_tmtpro_10855 | develop | [2405, 2407, 2443] | [2293, 2331, 2294] | False |
| eclipse_tmtpro_10855 | #10377 | [2405, 2407, 2443] | [2405, 2407, 2443] | True |

### ProSE TOPP tests (replayed through `ProSEAlgorithm::search()`)

The harness cannot link the TOPP tool (the wheel's libOpenMS has no TOPPBase). `topp_emulate.py` therefore replays
ProSE TOPP tests 1, 2, 4, 5, 6, 8 and 10 with their INI and command-line `Search:*` parameters. The `develop` replay
reproduces the checked-in `ProSE_1/2/4/10_out.idXML` SearchParameters and hits exactly.
- **#10377:** no change.
- **#10378:** one added search-parameter line.
- **#10379:** two lines (scorer) plus three lines (kernel). Hits and tolerances are unchanged.

The PRs' reference updates come from these replays (`update_refs.py`). CI confirmed them on Linux and macOS.

## Interpretation and limits

- **#10377** is yield-neutral on raw data and removes the intensity-scale dependence (scaled = unscaled, byte-identical PIN).
- **#10378:** the priors are opt-in.
  - Gains on Astral (+4.6%), timsTOF (+2.4%) and TMTpro (+1.8%).
  - About ±0.4% elsewhere; Velos CID −0.4%.
  - Native scores are unchanged.
  - The model trains in-file before Percolator cross-validation; this is not an out-of-fold assessment.
- **#10379, default change (`calibration:enabled=auto`):** this changes results for every high-resolution search.
  - Gains: Astral +5.1%, timsTOF +0.8%, TMTpro +0.6%.
  - Losses on unlabelled HCD: HF-X −3.3%, Lumos LFQ −4.2%.
  - Calibration tightens the precursor window on these files: `lumos_lfq_5192` goes from ±20 ppm to −3.50/+3.11 ppm,
    the same window as on 2026-09-30, and `hfx_A2` to ≤ 8.3 ppm. The 2026-09-30 diagnostic attributed most of the lost
    assignments to isotope-offset precursors outside that window; that diagnostic was not re-run here.
  - The opt-in scorer without calibration is neutral to positive (Astral +4.1%).
  - The fitted kernel itself is neutral (−1.2%…+0.7% vs. the fixed kernel) and lowers native TDC strongly on HCD.
- **#10335 merged with develop:**
  - CID gains remain: Velos +11.7%, Lumos CID TMT +3.0% (calibrated scorer with multiple fragment charges).
  - High-resolution groups are now within ±0.3%, because develop already has the window fix.
  - timsTOF −1.3%.
- **Scope:** counts are nominal-q yields, not validated true positives, and seeds are algorithmic repeats. The scope is
  sampled DDA from two vendors and six instrument models. This is not a full TOPP or pyOpenMS build; the pyOpenMS
  binding sources were only syntax-checked against nanobind 3.1.0.

## Reproduce

Place `pr-untangled/` (from `pr-untangled-reproduction.zip`) next to `prose-andes-reproduction/`. Prepare inputs, the
nightly and Percolator as in the original bundle README; Percolator also needs `libboost-filesystem1.83.0`.
Use pyOpenMS `3.7.0.dev20261001` (with `nanobind-backend` in `deps`). Check out the commits in `sources.json` as
worktrees `pr-untangled/sources/<build>`; for `pr10335`, apply `pr10335-on-develop.patch` to develop `f5ea2d04`. Then:

```bash
python3 pr-untangled/build.py base pr10368 pr10366 pr10365_scorer pr10365 pr10335
python3 pr-untangled/run.py --workers 2
python3 pr-untangled/intensity_control.py prepare
python3 pr-untangled/intensity_control.py run
python3 pr-untangled/intensity_control.py evaluate
python3 pr-untangled/evaluate.py --require-complete
python3 pr-untangled/topp_emulate.py base pr10365 --openms pr-untangled/sources/base   # TOPP replays
```

`reference/results/<dataset>/<arm>/` holds every job's parameters, resolved search parameters, input checks, command
records, search/export/Percolator logs and per-seed results with output hashes. Large PIN/TSV/idXML outputs and
binaries are regenerated.
