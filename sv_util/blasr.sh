#!/bin/bash

~/workspace/legacy_blasr/blasr $1 $2 -printSAMQV -nproc 32 -sam -out $3 > /dev/null