#!/bin/bash


mysnakefile="DIMet/tests/Snakefile_with_metabologram.smk"

cd $HOME # or any location above DIMet
cp -r DIMet/examples/toy*/ $HOME

snakemake -s $mysnakefile --cores 1 --config PRIMARY_CONFIG="~/toy_metabologram/analysis001/config-2-001.yml" \
METABOLOGRAM_CONFIG="~/toy_metabologram/analysis001/config_metabologram.yml" --latency-wait 15

