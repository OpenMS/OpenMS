#!/bin/sh
# The runs of this benchmark, in order. OPENMS_BIN: PR build (349328c1); PROTO: extracted prose-andes-reproduction bundle.
H=$(dirname "$0")
python3 -I $H/run.py std
python3 -I $H/run.py hardmods hfx_A2 hfx_B1 hfx_B3 astral_A2 astral_B1 astral_B3 eclipse_tmtpro_10855 eclipse_tmtpro_10858 eclipse_tmtpro_10863 velos_125_R1 velos_5000_R2 velos_25000_R3
python3 -I $H/run.py iso astral_A2 astral_B1 astral_B3 hfx_A2 hfx_B1 hfx_B3
python3 -I $H/run.py semi hfx_A2 hfx_B1 hfx_B3 velos_125_R1 velos_5000_R2 velos_25000_R3 astral_A2 astral_B1 astral_B3
python3 -I $H/run.py immuno liver_classI_r1
# follow-up: external Percolator (PERCOLATOR=<3.06 binary>) on the 20 files; isoforms on whole runs (fetch_full.sh first)
python3 -I $H/run.py ext
python3 -I $H/run.py isofull hfx_A2_full eclipse_tmtpro_10855_full astral_A2_full
python3 -I $H/evaluate.py ${WORK:-work} results
