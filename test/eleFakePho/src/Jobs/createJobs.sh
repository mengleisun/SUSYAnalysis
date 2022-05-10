#!/bin/bash

era=$1
nfile=$2
mkdir -p "condor/"$era

sed -e "s@ERA@$era@" run.sh.tmpl > "run_"$era".sh"
chmod 775 "run_"$era".sh"
sed -e "s@ERA@$era@" -e "s@NFILE@$nfile@" batch.jdl.tmpl > "batch_"$era".jdl"
