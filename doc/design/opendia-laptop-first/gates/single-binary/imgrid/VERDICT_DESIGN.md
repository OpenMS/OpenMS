# PLAN_IMGRID design verdict (committed before the held-out run)

Plan sha256 check: OK. NEW = 963543c7f0 (`OpenDIAWorkflow.imgrid`), BASE = 6dd893b954. `EVAL_D_img.txt`.

| rule | BASE | NEW | result |
|---|---|---|---|
| 1. D median N ≥ 1.03 × BASE, min NEW > max BASE | 14,875 (14,804-14,927) | 16,116 (15,957-16,177), +8.3 % | PASS |
| 2. D agreement ≥ BASE − 0.005 | 0.9344 | 0.9299 (−0.0045) | PASS (close to the limit) |
| 3. Astral raw01 byte-identical | | identical | PASS |

Design verdict: PASS. Held-out H decides rule 4.
