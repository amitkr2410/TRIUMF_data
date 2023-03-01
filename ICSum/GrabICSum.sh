#!/bin/bash
echo "Processing odb file (ICScalar Free Trigger Accepted Trigger)"
source /home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/ICScaler.sh
/bin/bash /home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/ICScaler.sh

echo "Processing ICScalar.C"
root -b -q "/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/ICScalar.C"

echo "Processong ICScalarPlot.C"
root -b -q "/home/iris/Amit/Analysis_10C_8torr/Plot/ICSum/ICScalarPlot.C"
