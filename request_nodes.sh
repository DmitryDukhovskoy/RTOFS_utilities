#!/bin/bash
# request comput node for interactive session
# to run python or debugging
# max hours - 12 but check:
# sacctmgr show qos format=Name,MaxWall 
#salloc --x11 -q batch -t 10:00:00 --ntasks=1 --clusters=c5 -A cefi
# On gaea:
salloc -M c5 -t 10:00:00
