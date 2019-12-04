#!/bin/bash
set -e # exit on error
while [ 1 ]
do
    python3 create_sample_dataset.py
    python3 compare_callers.py
    rm -r /MAdata/sv_datasets/minimal-3/
done