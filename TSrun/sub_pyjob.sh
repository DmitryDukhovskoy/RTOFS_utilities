#!/bin/sh 
#
# Submit python job 
# on interactive node
# sbatch -q batch -t 8:00:00 --nodes=1 -A marine-cpu sub_pyjob.sh
#
export WD=/scratch2/NCEPDEV/marine/Dmitry.Dukhovskoy/python/TSrun
cndn=anls
pexe=tsprof_error_anls.py

cd $WD
pwd
ls -la
eval "$($PYPATH/bin/conda shell.bash hook)"
conda activate $cndn

which ipython

echo "Submitting ipython < ${pexe}"
ls -l ${pexe}
#ipython < run ${pexe}
ipython ${pexe} >> ${pexe}.log

exit 0

