#!/bin/bash

echo "This is UnlabeledMatcher.sh - a small shell script to illustrate the usage of UnlabeledMatcher."
echo "The feature map data files can be obtained upon request from the authors."

#--------------------------------------------------
function cmd-wrapper()
{
# uncomment your choice (for debugging of this script)
#  yes, DO run
	time $@
#  simulate
#	echo cmd-wrapper: $@
}
#--------------------------------------------------


echo " "
echo "============================================================"
echo " "

echo "We start with a simple, artificial example.  Two maps, shifted by 3 seconds in retention time."
echo "As was to be expected, both algorithms find all 2521 pairs."
echo " "

cmd-wrapper UnlabeledMatcher -ini UnlabeledMatcher.ini -n 1 -d 1000
cmd-wrapper UnlabeledMatcher -ini UnlabeledMatcher.ini -n 2 -d 1000
# time UnlabeledMatcher -ini UnlabeledMatcher.ini -n 3 -d 1000

echo " "
echo "============================================================"
echo " "

echo "The next example shows two repeat measurements."
echo "The 'geomhash_shift' algorithm finds 1928, the 'simple' algorithm 1890 feature pairs."
echo " "

cmd-wrapper UnlabeledMatcher -ini UnlabeledMatcher.ini -n 10 -d 1
cmd-wrapper UnlabeledMatcher -ini UnlabeledMatcher.ini -n 11 -d 1

cmd-wrapper gnuplot E32R_24h_1_vs_2.geomhash_shift.data.gp - <<EOF
set output 'E32R_24h_1_vs_2.geomhash_shift.data.png'
set terminal png size 1000,1000
replot
quit
EOF

cmd-wrapper gnuplot E32R_24h_1_vs_2.simple.data.gp - <<EOF
set output 'E32R_24h_1_vs_2.simple.data.png'
set terminal png size 1000,1000
replot
quit
EOF

echo " "
echo "The difference in the number of feature pairs is not so spectacular,"
echo "however the qualities of the feature pairs are much better,"
echo "as is illustrated by the plot we generate next."
echo "The x-axis is the index of the feature pair, ordered by quality,"
echo "the y-axis is the quality."
echo " "

cmd-wrapper FeaturePairSplitter -in E32R_24h_1_vs_2.geomhash_shift.pairs.xml -qual E32R_24h_1_vs_2.geomhash_shift.pairs.qualities
LC_ALL=C sort -nr E32R_24h_1_vs_2.geomhash_shift.pairs.qualities > E32R_24h_1_vs_2.geomhash_shift.pairs.qualities.sorted

cmd-wrapper FeaturePairSplitter -in E32R_24h_1_vs_2.simple.pairs.xml -qual E32R_24h_1_vs_2.simple.pairs.qualities
LC_ALL=C sort -nr E32R_24h_1_vs_2.simple.pairs.qualities > E32R_24h_1_vs_2.simple.pairs.qualities.sorted

cmd-wrapper gnuplot - <<EOF
plot 'E32R_24h_1_vs_2.geomhash_shift.pairs.qualities.sorted' with lines
replot 'E32R_24h_1_vs_2.simple.pairs.qualities.sorted' with lines
set output 'E32R_24h_1_vs_2.algorithm.qualities.png'
set term png
replot
quit
EOF

echo " "
echo "============================================================"
echo " "

echo "Finally, we go through the same steps for another two repeat measurements."

echo "Here the 'geomhash_shift' algorithm finds 2182, the 'simple' algorithm 2155 feature pairs."
echo "Have a look at the  grid file.  It lists the shifts,"
echo "which change over time in a complex, non-linear way."
echo "The simple algorithm creates quite a few false pairs in the cluster near RT=500"
echo " "

cmd-wrapper UnlabeledMatcher -ini UnlabeledMatcher.ini -n 20 -d 1
cmd-wrapper UnlabeledMatcher -ini UnlabeledMatcher.ini -n 21 -d 1

cmd-wrapper gnuplot E32R_24h_10uMCu_1_vs_2.geomhash_shift.data.gp - <<EOF
set output 'E32R_24h_10uMCu_1_vs_2.geomhash_shift.data.png'
set terminal png size 1000,1000
replot
quit
EOF

cmd-wrapper gnuplot E32R_24h_10uMCu_1_vs_2.simple.data.gp - <<EOF
set output 'E32R_24h_10uMCu_1_vs_2.simple.data.png'
set terminal png size 1000,1000
replot
quit
EOF

cmd-wrapper FeaturePairSplitter -in E32R_24h_10uMCu_1_vs_2.geomhash_shift.pairs.xml -qual E32R_24h_10uMCu_1_vs_2.geomhash_shift.pairs.qualities
LC_ALL=C sort -nr E32R_24h_10uMCu_1_vs_2.geomhash_shift.pairs.qualities > E32R_24h_10uMCu_1_vs_2.geomhash_shift.pairs.qualities.sorted

cmd-wrapper FeaturePairSplitter -in E32R_24h_10uMCu_1_vs_2.simple.pairs.xml -qual E32R_24h_10uMCu_1_vs_2.simple.pairs.qualities
LC_ALL=C sort -nr E32R_24h_10uMCu_1_vs_2.simple.pairs.qualities > E32R_24h_10uMCu_1_vs_2.simple.pairs.qualities.sorted

cmd-wrapper gnuplot - <<EOF
plot 'E32R_24h_10uMCu_1_vs_2.geomhash_shift.pairs.qualities.sorted' with lines
replot 'E32R_24h_10uMCu_1_vs_2.simple.pairs.qualities.sorted' with lines
set output 'E32R_24h_10uMCu_1_vs_2.algorithm.qualities.png'
set term png
replot
quit
EOF

echo " "
echo "============================================================"
echo " "

echo " "
echo "Done."
echo " "


