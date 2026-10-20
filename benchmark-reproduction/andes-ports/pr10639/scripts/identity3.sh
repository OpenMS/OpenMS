#!/bin/bash
# Slices off vs. 4 / 16 with the total peptide order (/opt/openms-v3): sha256 of every PIN and of every idXML
# without its date / database path lines. One line per case and arm.
P=/opt/openms-v3; D=/home/user/bench/frozen/z/prose-andes-reproduction/broad_benchmark/data; R=runs/v3_identity
h() { for f in "$1"/*.pin; do sha256sum "$f" | cut -c1-16; done | tr '\n' ' '
      for f in "$1"/*.idXML; do grep -v 'date=\|db=\|file_origin' "$f" | sha256sum | cut -c1-16; done | tr '\n' ' '; rm -f "$1"/*.idXML; }
run() { case=$1 ini=$2 db=$3; shift 3; in=("$@")
  for s in $SL; do o=$R/$case/s$s; scripts/run1.sh $P $ini $o 4 $db "${in[@]}" -- -Search:database:mass_slices $s >/dev/null && echo "$case s$s $(h $o)"; done; }
SL="1 4 16" run frozen ini/hfx_A2.ini db/hye_E.fasta $D/hfx_A2.mzML
SL="1 4 16" run multi3 ini/hfx_A2.ini db/hye_E.fasta $D/hfx_A2.mzML $D/hfx_B1.mzML $D/astral_B1.mzML
SL="1 4 16" run whole ini/hfx_A2.ini db/hye_E.fasta whole/hfx_A2_full.mzML
SL="4 16" run phos2_E ini/hfx_A2_phos2.ini db/hye_E.fasta $D/hfx_A2.mzML
SL="4 16" run mhc ini/mhc.ini mhc/human_td.fasta mhc/DN08_Brain_classI_techRep1.mzML
echo IDENTITY_DONE
