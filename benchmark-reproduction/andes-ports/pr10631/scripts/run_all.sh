#!/bin/bash
# Reproduces the PR #10631 benchmark. Needs three OpenMS installs (bin/ lib/ of a WITH_THERMO_RAW=ON build) in
# ${OPENMS_PREFIX:-/opt}/openms-{base,pr,fix}, the inputs in $DATA (see ../README.md) and python3 with numpy.
set -e
S=$(cd "$(dirname "$0")" && pwd); export DATA=${DATA:-/home/user/data}; cd $DATA
mkdir -p search quant
# 1. .raw -> mzML (the converter is not touched by the PR; one build for all arms)
source $S/env.sh pr
for s in Fig2_m23141_hrMS2 Fig2_m23142_hrMS2 Fig2_m23145_SPSms3 Fig2_m23146_SPSms3; do
  FileConverter -in pxd017803/$s.raw -out pxd017803/$s.mzML -threads 2
done
for f in pxd017803/*.mzML arch/*.mzML; do python3 -I $S/scan_table.py $f ${f%.mzML}.scans.tsv; done
# 2. identification, once (identical in all arms): ProSE, in-process Percolator, 1% PSM FDR
$S/search_tko.sh Fig2_m23145_SPSms3 0.5 Da 2; $S/search_tko.sh Fig2_m23146_SPSms3 0.5 Da 2
$S/search_tko.sh Fig2_m23142_hrMS2 20 ppm 2;  $S/search_tko.sh Fig2_m23141_hrMS2 20 ppm 2
$S/search_arch.sh a05058 0.5 Da TMT6plex 2;   $S/search_arch.sh xb10855 20 ppm TMTpro 2
for f in search/*.idXML; do python3 -I $S/idxml_table.py $f ${f%.idXML}.ids.tsv; done
# 3. quantification per arm: IsobaricAnalyzer and IsobaricWorkflow, default parameters
declare -A TYPE=([a05058]=tmt11plex [xb10855]=tmt18plex)
for arm in base pr fix; do
  for s in Fig2_m23145_SPSms3 Fig2_m23146_SPSms3 Fig2_m23142_hrMS2 Fig2_m23141_hrMS2; do
    $S/run_quant.sh $arm pxd017803/$s.mzML search/$s.idXML tmt16plex quant/$arm; done
  for s in a05058 xb10855; do $S/run_quant.sh $arm arch/$s.mzML search/$s.idXML ${TYPE[$s]} quant/$arm; done
  # single survey scan (no interpolation), for the comparison with the reference purity
  mkdir -p quant/${arm}_noint; source $S/env.sh $arm
  for s in a05058 xb10855 Fig2_m23145_SPSms3 Fig2_m23142_hrMS2; do
    sub=arch; t=${TYPE[$s]:-tmt16plex}; [ -f pxd017803/$s.mzML ] && sub=pxd017803
    IsobaricWorkflow -in $sub/$s.mzML -in_id search/$s.idXML -out quant/${arm}_noint/$s.IW.consensusXML -type $t \
      -calculate_id_purity -extraction:purity_interpolation false > quant/${arm}_noint/$s.IW.log 2>&1
  done
done
for f in quant/*/*.consensusXML; do python3 -I $S/cons_table.py $f ${f%.consensusXML}.tsv; done
# 4. evaluation
for s in a05058 xb10855; do python3 -I $S/ref_purity.py arch/$s.mzML arch/$s.scans.tsv arch/$s.refpurity.tsv; done
for s in Fig2_m23145_SPSms3 Fig2_m23142_hrMS2; do python3 -I $S/ref_purity.py pxd017803/$s.mzML pxd017803/$s.scans.tsv pxd017803/$s.refpurity.tsv; done
python3 -I $S/evaluate.py $DATA quant_summary.tsv Fig2_m23145_SPSms3:tko@TKO_SPS-MS3 Fig2_m23146_SPSms3:tko@TKO_SPS-MS3 \
  Fig2_m23142_hrMS2:tko@TKO_hrMS2 Fig2_m23141_hrMS2:tko@TKO_hrMS2 a05058 xb10855
python3 -I $S/compare_ref.py $DATA purity_vs_reference.tsv a05058 xb10855 Fig2_m23145_SPSms3 Fig2_m23142_hrMS2
python3 -I $S/spikein_ratios.py $DATA spikein_ratios_a05058.tsv
for f in pxd017803/*.scans.tsv arch/*.scans.tsv; do python3 -I $S/scan_stats.py $f; done > acquisition.txt
