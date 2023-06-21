This example folder provides examples of datasets:

Copy and them in your 'home/'

toy1 : time course example
toy2 : two conditions comparison
toy3 : two conditions comparison with special grouping
toy_multigroup : >= 3 conditions being compared 

For example, for running toy1, in your home, type:

```
conda activate dimet
python -m DIMet.src.prepare $HOME/toy1/analysis001/config-1-001.yml
```

The 'toy_metabologram' reuses the toy2 data, but uses 
selected results for running the metabologram.

