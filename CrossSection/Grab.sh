#!/bin/bash
for ((  i = 1;  i < 17;  i++  ))
do
 echo "...........Analyzing Ring = $i ......"
 root -b -q "QValueEachRing.C($i)" 

done
