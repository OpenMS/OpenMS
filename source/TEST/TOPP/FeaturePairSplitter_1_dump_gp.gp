# (all comments removed)
# (all comments removed)
plot   "FeaturePairSplitter_1_dump.tmp" using 2:3 title "map 1"
replot "FeaturePairSplitter_1_dump.tmp" using 5:6 title "map 2"
replot "FeaturePairSplitter_1_dump.tmp" using 2:3:($5-$2):($6-$3) w vectors nohead title "pairs"
# (all comments removed)
