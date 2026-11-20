# RESULT_TRPPM_H: held-out run reproduces the gain (+9.9 % IDs at 1 % true FDP). Reported after the D verdict

- **Order:** D verdict committed as 01555b47df at 15:50:42 UTC, and H submitted at 15:50:46 (jobs 38728635-44).
- **Binary and setup:** as in PLAN_TRPPM (frozen 0d80832844).
- **Evaluation:** `eval_pdec.py` (bb03387d…). Output is in `trppm/EVAL_H.txt`, and wall time and RSS are in `trppm/fetch_H.txt`.

| Arm | Seed | N @ q ≤ 1 % | E | Lower | Combined | N @ FDP_c 1 % |
|---|---|---|---|---|---|---|
| CUR | 0 | 27,274 | 300 | 0.0110 | 0.0126 | 26,697 |
| CUR | 1 | 27,280 | 289 | 0.0106 | 0.0121 | 26,577 |
| CUR | 2 | 27,308 | 292 | 0.0107 | 0.0122 | 26,647 |
| CUR | 3 | 27,071 | 273 | 0.0101 | 0.0115 | 26,651 |
| CUR | 4 | 27,180 | 287 | 0.0106 | 0.0121 | 26,522 |
| NEW | 0 | 30,003 | 367 | 0.0122 | 0.0140 | 29,236 |
| NEW | 1 | 30,119 | 346 | 0.0115 | 0.0131 | 29,286 |
| NEW | 2 | 30,084 | 356 | 0.0118 | 0.0135 | 29,434 |
| NEW | 3 | 29,996 | 340 | 0.0113 | 0.0129 | 29,353 |
| NEW | 4 | 29,586 | 316 | 0.0107 | 0.0122 | 29,209 |

| Pooled | CUR | NEW |
|---|---|---|
| median N @ FDP_c 1 % | 26,647 | **29,286** (+9.9 %) |
| lower bound at own q ≤ 1 % | 0.0106 (90 % 0.0101-0.0111) | 0.0115 (0.0111-0.0120) |
| combined bound | 0.0121 | 0.0131 |
| pair win | 0.492 | 0.501 |
| median wall | 11:24 | 11:11 |
| peak RSS | ≤ 28.6 GB | ≤ 27.8 GB |

- **The gain holds on the held-out run.** Every NEW seed is above every CUR seed. The gain is smaller than on D (+9.9 % against +18.4 %), but it is still the largest single step on H in this series.
- **Calibration is a little looser on H.** The lower bound rises +0.0009 and the combined bound +0.0010, which is within what the D criteria allowed (+0.0030 / +0.0040). The tool's q ≤ 1 % on H is now about 1.15 % true FDP (lower bound). That is still part of the residual excess (b9).
- **OpenDIA against DIA-NN at FDP_c 1 %:** 29,286 against 41,272 on H, which is 71 % (it was 65 %).
