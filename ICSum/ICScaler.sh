#!/bin/bash

THISDIR=$PWD
rm "/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar.txt"
FILE="/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/Scalar.txt"
echo "#run ICScalar TotalFreeTriger TotalAcceptedTriger " >> $FILE
cd /ladd/iris_data0/iris/oct2013_S1396/

##/bin/cat <<EOF >$FILE

for (( i = 2178 ; i <2328 ; i++ )) do
 ICScalarTotal=`odbhist -r $i $i -v "/Equipment/Adcscaler/variables/scas[30]" -q`  ##For ICScalar = scas[30]
 FreeTrigTotal=`odbhist -r $i $i -v "/Equipment/Adcscaler/variables/scas[31]" -q`  ##For ICScalar = scas[31]
 AcceptTrigTotal=`odbhist -r $i $i -v "/Equipment/Adcscaler/variables/scas[26]" -q`  ##For ICScalar = scas[26]
 echo "Processing run number " $i  
 echo "No of beam particle" $a
 printf "\n"
echo $i $ICScalarTotal $FreeTrigTotal $AcceptTrigTotal >> $FILE
done
cd $THISDIR
