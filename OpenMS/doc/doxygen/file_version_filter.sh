#!/bin/sh
grep -o -e "\$Maintainer: .*\$" $1 | sed 's/\$//g'
