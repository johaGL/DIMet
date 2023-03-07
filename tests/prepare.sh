#!/bin/bash

# the test needs toys to be in $HOME
echo $PWD

cd $HOME

echo $PWD

echo "===> toy1"

python -m DIMet.src.prepare $HOME/toy1/analysis001/config-1-001.yml


echo "===> toy2"
echo "config 1"
python -m DIMet.src.prepare $HOME/toy2/analysis001/config-2-001.yml
echo "config 2"
python -m DIMet.src.prepare $HOME/toy2/analysis002/config-2-002.yml --use_internal_standard Myristic_acid_d27


echo "===> toy3"

python -m DIMet.src.prepare $HOME/toy3/analysis01/config-3-01.yml



