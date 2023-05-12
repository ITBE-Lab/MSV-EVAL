#!/bin/bash
~/workspace/bwa/bwa mem -R \"@RG\\tID:1\\tSM:$1\" -t 32 $2 $3 > $4 2> /dev/null