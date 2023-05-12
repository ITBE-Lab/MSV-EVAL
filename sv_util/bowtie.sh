#!/bin/bash
~/workspace/bowtie2/bowtie2-2.3.3.1/bowtie2 --rg-id 1 --rg SM:$1 -p 32 -x $2 -1 $3 -2 $4 -S $5 2> /dev/null