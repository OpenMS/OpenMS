**A. #10379 windows, HyperScore. % vs develop / native TDC delta vs develop**

| Arm | HF-X HCD | Astral HCD | Lumos HCD LFQ | Exploris 480 TMTpro | timsTOF HT |
| --- | --- | --- | --- | --- | --- |
| #10379 head: `auto` = precursor + fragment windows | -3.34% / -4 | +5.09% / +331 | -4.16% / -142 | +0.61% / +174 | +0.79% / -4 |
| precursor window only | -2.83% / -9 | -3.46% / -5 | -3.71% / -148 | -0.87% / -2 | +0.79% / -4 |
| fragment window only | -0.91% / +27 | +5.48% / +349 | -0.43% / +15 | +1.94% / +206 | +0.00% / +0 |
| robust precursor window only (ANDES-style) | -4.73% / -68 | -5.15% / +159 | -7.24% / -302 | -5.97% / -67 | -0.57% / -6 |
| robust precursor window, stride sample | -6.40% / -119 | -4.44% / +167 | -7.09% / -287 | -7.34% / -69 | -3.12% / -26 |

| (develop PSMs) | 4380.3 | 2189.6 | 5496.1 | 2244.9 | 1015.0 |

**B. Mass-accuracy scorer. % vs develop (vs scorer with fixed 7 ppm kernel, no calibration) / native TDC delta vs develop**

| Arm | HF-X HCD | Astral HCD | Lumos HCD LFQ | Exploris 480 TMTpro | timsTOF HT |
| --- | --- | --- | --- | --- | --- |
| scorer, fixed 7 ppm kernel, no calibration | +0.34% (+0.00) / -13 | +4.13% (+0.00) / +417 | -0.14% (+0.00) / +37 | -0.07% (+0.00) / +231 | +0.62% (+0.00) / +10 |
| #10379 head: fitted kernel + both windows | -3.04% (-3.37) / -601 | +3.10% (-0.98) / -186 | -3.39% (-3.26) / -269 | +1.22% (+1.29) / +72 | +1.25% (+0.62) / +48 |
| fitted kernel (center + width), windows unchanged | +0.25% (-0.09) / -714 | +0.42% (-3.56) / -204 | -0.17% (-0.03) / -194 | -0.07% (-0.00) / +52 | -0.72% (-1.34) / +52 |
| fitted center only, windows unchanged | +0.09% (-0.25) / +10 | +4.95% (+0.79) / +383 | -0.08% (+0.06) / +35 | -0.45% (-0.38) / +227 | -0.57% (-1.19) / +61 |
| fitted center only + fragment window | -0.85% (-1.19) / +12 | +6.52% (+2.30) / +404 | -0.44% (-0.30) / +18 | +2.22% (+2.29) / +232 | -0.57% (-1.19) / +61 |
| fitted center only + robust precursor window, stride | -5.39% (-5.71) / -138 | -0.97% (-4.90) / +412 | -7.15% (-7.02) / -312 | -8.10% (-8.03) / +86 | -3.96% (-4.55) / +22 |

**C. #10378 ion priors. % vs develop (vs #10378 head)**

| Arm | Velos CID | HF-X HCD | Astral HCD | Lumos HCD LFQ | Lumos CID TMT | Exploris 480 TMTpro | timsTOF HT |
| --- | --- | --- | --- | --- | --- | --- | --- |
| #10378 head (fragment charge 1) | -0.37% (+0.00) | +0.30% (+0.00) | +4.60% (+0.00) | -0.11% (+0.00) | +0.92% (+0.00) | +1.76% (+0.00) | +2.35% (+0.00) |
| fragment charges `auto` | +6.20% (+6.59) | +0.31% (+0.00) | – | – | +3.14% (+2.19) | – | – |
| cleavage-residue context | +0.54% (+0.91) | +0.63% (+0.33) | +4.82% (+0.20) | -0.04% (+0.07) | +1.92% (+0.99) | +1.44% (-0.32) | +1.72% (-0.61) |
| both | +6.16% (+6.55) | – | – | – | +4.03% (+3.08) | – | – |

