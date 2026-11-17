# PLAN_RTCENTRE design verdict (committed before the held-out run)

Plan sha256 check: OK. Runs: `run_rtc.sh`, evaluation `eval_rtc.py` -> `EVAL_D.txt`.

| rule | BASE (0e9e4fc3a5) | NEW (6dd893b954) | result |
|---|---|---|---|
| 1. D median N ≥ 1.02 x BASE, min NEW > max BASE | median 13,127 (13,109-13,290) | median 14,875 (14,804-14,927), +13.3 % | PASS |
| 2. D agreement with DIA-NN ≥ BASE − 0.005 | 0.9308 | 0.9344 | PASS |
| 3. D map p95 NEW < BASE | 19.31 s (early median −19.4 s) | 4.73 s (early +0.1 s) | PASS |
| 4. Astral raw01 byte-identical | | precursors.tsv, run precursors.tsv, calibration.json identical | PASS |

Design verdict: PASS. Held-out H (HeLa 22474.d) is run next and decides rule 5.
