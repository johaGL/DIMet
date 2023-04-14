#!/bin/bash

cd $HOME # or any location above DIMet

snakefi="DIMet/tests/Snakefile_yes_isotopolo.smk"

snakemake -s $snakefi --cores 1 --config PRIMARY_CONFIG="~/toy1/analysis001/config-1-001.yml"
#METABOLOGRAM_CONFIG="~/toy1/analysis001/metabologram_config.yml" --latency-wait 15

echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""

snakemake -s $snakefi --cores 1 --config PRIMARY_CONFIG="~/toy2/analysis001/config-2-001.yml"
echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""

snakemake -s $snakefi --cores 1 --config PRIMARY_CONFIG="~/toy2/analysis002/config-2-002.yml"
echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""

snakemake -s $snakefi --cores 1 --config PRIMARY_CONFIG="~/toy3/analysis01/config-3-01.yml"
echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""
