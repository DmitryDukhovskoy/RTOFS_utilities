#!/bin/sh -x
# Swap dir names in python scripts
# 
set -u

#DR=/home/Dmitry.Dukhovskoy/python/rtofs
DR=/home/Dmitry.Dukhovskoy/python/prepare_relax
DO=/scratch2/NCEPDEV/marine
DN=/home
cd $DR
pwd

for fl in *.py
do
  echo $fl
  sed -i "s|^sys.path.append('${DO}|sys.path.append('${DN}|g" $fl
done

