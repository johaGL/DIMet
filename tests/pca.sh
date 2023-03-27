#!/bin/bash

# the test needs toys to be in $HOME
echo $PWD

cd $HOME

echo $PWD

echo "===> toy1"
python -m DIMet.src.pca $HOME/toy1/analysis001/config-1-001.yml --run_iris_demo


echo "===> toy2"
echo "config 1"
python -m DIMet.src.pca $HOME/toy2/analysis001/config-2-001.yml
echo "config 2"
python -m DIMet.src.pca $HOME/toy2/analysis002/config-2-002.yml


echo "===> toy3"

python -m DIMet.src.pca $HOME/toy3/analysis01/config-3-01.yml



