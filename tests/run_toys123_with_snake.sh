#!/bin/bash

snakef_common="DIMet/tests/Snakefile_common.smk"

cd $HOME 
cp -r DIMet/examples/toy*/ $HOME

snakemake -s $snakef_common --cores 1 --config PRIMARY_CONFIG="~/toy1/analysis001/config-1-001.yml"


echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""

snakemake -s $snakef_common --cores 1 --config PRIMARY_CONFIG="~/toy2/analysis001/config-2-001.yml"
echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""

snakemake -s $snakef_common --cores 1 --config PRIMARY_CONFIG="~/toy2/analysis001/config-2-001.yml"
echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""

snakemake -s $snakef_common --cores 1 --config PRIMARY_CONFIG="~/toy3/analysis001/config-3-01.yml"
echo ""
echo " ** ** ** ** ** ** ** ** ** ** **"
echo ""
