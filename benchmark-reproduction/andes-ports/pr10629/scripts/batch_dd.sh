# DecoyDatabase -entrapment_ratio arms (PR build): $3 = arm suffix, $4 = hye database, $5 = human database
S=$1; B=$2; ARM=$3; HYE=$4; HUMAN=$5
cd $S
for d in hfx_A2 hfx_B1 hfx_B3 astral_A2 astral_B1 astral_B3; do python3 -I scripts/run.py $B $d $HYE runs/$d/$ARM; done
[ -n "$HUMAN" ] && for d in eclipse_tmtpro_10855 eclipse_tmtpro_10858 eclipse_tmtpro_10863; do python3 -I scripts/run.py $B $d $HUMAN runs/$d/$ARM; done
echo BATCH_DONE $ARM
