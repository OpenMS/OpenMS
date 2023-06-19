#!/bin/bash
# Use in your build folder with the test regex (as for ctest -R) as argument like: ../tools/overwriteTOPPTestFiles.sh TOPP_OpenSwathAnalyzer
ctest -V -R $1 > tmp.tst.log
grep -e "^[0-9]*: diff" tmp.tst.log | sed 's/^.*diff/cp/g' > copies_to_perform.txt
cat copies_to_perform.txt | while read -r line
do
   echo "$line"
   $line
done
rm tmp.tst.log copies_to_perform.txt