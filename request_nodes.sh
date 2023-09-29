#!/bin/bash
# request comput node for interactive session
# to run python or debugging
# max hours - 8 (?)
salloc --x11=first -q batch -t 8:00:00 --nodes=1 -A marine-cpu
