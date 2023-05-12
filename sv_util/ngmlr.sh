#!/bin/bash
$1 -r $2 -q $3 --rg-id 1, --rg-sm $4 -t 8 -x $5 > $6 || rm $6 
#2>/dev/null
