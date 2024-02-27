#!/bin/bash -vx
## Need to request compute node for interactive session first
## to run python or debugging
## max hours - 8 (?)
## start conda and activate session
cndn="anls"

# Request compute node:
#./request_nodes.sh

#echo "status =$?"
#if [[ $?=0 ]]; then 
  eval "$($PYPATH/bin/conda shell.bash hook)"
  conda activate $cndn
  which ipython
  echo "python $cndn is active, ready to run"
#else
#  echo "Could not get compute node"
#  exit 1
#fi

