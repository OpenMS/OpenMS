#!/bin/bash


for i in `seq 1 31`;
do
	UnlabeledMatcher -ini AddSeries.ini -n $i
	MapMatcher -ini AddSeries.ini -n $i
	Dewarper -ini AddSeries.ini -n $i
		
done  

AdditiveSeries -ini AddSeries.ini 
    
