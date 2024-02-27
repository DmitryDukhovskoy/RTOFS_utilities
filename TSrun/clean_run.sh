#!/bin/sh -x
set -u

export WD=/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/TSrun
cd $WD
pwd
ls -l

# Clean temporary files:
/bin/rm -f derive_run*py
/bin/rm -f *.log
/bin/rm -f pyjob*.sh
/bin/rm -f slurm-*out
/bin/rm -f testrun*log

exit 0
