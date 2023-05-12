#!/bin/bash

# -c output CIGAR in PAF; -a output SAM
$1 --MD -c -a -t 32 -x $2 $3 -R @RG\\tID:1\\tSM:$4 $5 $6 > $7 2>/dev/null