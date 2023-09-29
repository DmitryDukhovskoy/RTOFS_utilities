#!/bin/sh -x
#
# Get QC text log files with T/S profiles
# Use for deriving Argo T/S profiles in text files
#
set -u


export DRUN=NCEPDEV
export expt=paraB

#/NCEPDEV/emc-ocean/5year/Dan.Iredell/wcoss2.prod/rtofs.[YYYYMMDD]
#first day is 20220526

obs="profile_qc"  # what obs to keep

if [[ $#<1 ]]; then
  printf " ERR: Usage get_qctar.sh YYYYMMDD [e.g., 20230123]"
  exit 1
fi
 
RD=$1
export D='/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/data'
export DUMP="${D}/rtofs.$RD"
mkdir -pv $DUMP
cd ${DUMP}

print "Checking HPSS"
hsi -P ls -l /${DRUN}/emc-ocean/5year/Dan.Iredell/wcoss2.${expt}/rtofs.${RD}/
print "Untarring logs tar"
#htar -xvf /NCEPPROD/5year/hpssprod/runhistory/rh${RD:0:4}/${RD:0:6}/$RD/com_rtofs_prod_rtofs.$RD.ab.tar ./'*'n-24.archv.'*' ./'*'n00.archv.'*'
htar -xvf /${DRUN}/emc-ocean/5year/Dan.Iredell/wcoss2.${expt}/rtofs.${RD}/logs.tar  

wait 

pwd
ls -l

/bin/rm -r ocnqc_logs
/bin/mv logs ocnqc_logs

cd ocnqc_logs
pwd
ls -l

for dir in */
do
  if [ ! $dir == "$obs/" ];  then
    printf " removing $dir "
    /bin/rm -r $dir
  fi
done 

# Remove everything but argo profiles:
if [ -d profile_qc ]; then
  cd profile_qc
  for fl in $(ls --ignore=*argo*)
  do
    /bin/rm $fl
  done
fi


exit 0

