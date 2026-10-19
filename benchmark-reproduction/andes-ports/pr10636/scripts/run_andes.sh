#!/bin/sh
# andes reference (bigbio/andes fa0195a, current main with #164) on the same frozen spectra: UniProt human reference
# proteome XML (UP000005640, release 2026_03) with --uniprot-isoforms, Percolator 3.06 (OpenMS THIRDPARTY) via --rescore.
# usage: run_andes.sh <andes binary> <frozen data dir> <UP000005640_9606.xml.gz> <out dir>
A=$1; D=$2; X=$3; O=$4; mkdir -p $O; cd $O
printf 'NumMods=1\n57.021463721,C,fix,any,Carbamidomethyl (C)\n15.994914622,M,opt,any,Oxidation (M)\n' > mods_none.txt
printf 'NumMods=1\n57.021463721,C,fix,any,Carbamidomethyl (C)\n15.994914622,M,opt,any,Oxidation (M)\n304.207146,K,fix,any,TMTpro (K)\n304.207146,*,fix,N-term,TMTpro (N-term)\n' > mods_TMTpro.txt
common="--database $X --uniprot-isoforms --threads 2 --charge 2..5 --top-n 10 --max-missed-cleavages 2 --min-length 7 --max-length 40 --max-mods 1 --precursor-cal off --isotope-error=-1..2 --rescore"
for f in astral_A2 astral_B1 astral_B3; do
  $A --spectrum $D/$f.mzML --output-pin $f.pin --output-tsv $f.tsv --precursor-tol 20ppm --mods mods_none.txt --protocol standard $common > $f.log 2>&1
done
for f in eclipse_tmtpro_10855 eclipse_tmtpro_10858 eclipse_tmtpro_10863; do
  $A --spectrum $D/$f.mzML --output-pin $f.pin --output-tsv $f.tsv --precursor-tol 50ppm --mods mods_TMTpro.txt --protocol TMT $common > $f.log 2>&1
done
