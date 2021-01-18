#!/bin/bash
ctest -V -R $1 > tmp.tst.log
grep -e "^[0-9]*: diff" tmp.tst.log | sed 's/^.*diff/cp/g' > copies_to_perform.txt
cat copies_to_perform.txt | while read line 
do
   echo "$line"
   $line
done
rm tmp.tst.log copies_to_perform.txt