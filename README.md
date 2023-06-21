# DIMet

[![License](https://img.shields.io/github/license/johaGL/dimet?label=license)](https://github.com/johaGL/dimet/blob/main/LICENSE)

[![CBiB Logo](imgs/cbib_logo.png)](https://www.cbib.u-bordeaux.fr/)

## DIMet: Differential Isotope-labeled targeted Metabolomics
----------------------------
**DIMet** is a bioinformatics pipeline for differential analysis of targeted isotope-labelled data.

DIMet supports the analysis of full metabolite abundances and isotopologue contributions, and allows to perform it either in the differential comparison mode or as a time-series analysis. As input, the DIMet accepts three types of measures: a) isotopologues’ contributions, b) fractional contributions (also known as mean enrichment), c) full metabolites’ abundances. Specific functions process each of the three types of measures separately.

**Note:** DIMet is intended for downstream analysis of tracer metabolomics data that has been **corrected** for the presence of natural isotopologues. Make sure you that the metabolomics platform provides you the output of the correction procedure before using this pipeline. 

![schema](imgs/schemaAlone4github.png)


## Requirements

You need a UNIX system, with conda or miniconda3 installed, see [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Then clone DIMet repository in your `$HOME` directory, and set the virtual environment by running (from terminal):
```
cd $HOME
git clone git@github.com:cbib/DIMet.git
conda env create --file DIMet/dimet.yml
conda activate dimet 
```

Before DIMet users may be interested in normalizations and other processing
that can be done with our acompanying tool https://github.com/johaGL/Tracegroomer :
- generate all the other tables using the isotopologue Absolute values table,  and/or
- normalize by the amount of material (number of cells, tissue weight), and/or
- normalize by an internal standard (present in your data) at choice

The output of Tracegroomer can be directly copied into your project to be analyzed by DIMet. 

# How to run DIMet

## Option 1 : With Jupyter notebooks

See [examples/jupyter_ldh.ipynb](examples/jupyter_ldh.ipynb) and
[jupyter_cyclos-time-course.ipynb](jupyter_cyclos-time-course.ipynb)


## Option 2 : With scripts in the unix terminal

Here we show you the raw scripts to analyze the data in
[examples/toy_metabologram/](examples/toy_metabologram/). 
Copy the entire example folder in your `$HOME`, and run:

```
python3 -m DIMet.src.prepare toy_metabologram/analysis001/config-2-001.yml
python3 -m DIMet.src.pca toy_metabologram/analysis001/config-2-001.yml
python3 -m DIMet.src.abundances_bars toy_metabologram/analysis001/config-2-001.yml
python3 -m DIMet.src.MEorFC_lineplot toy_metabologram/analysis001/config-2-001.yml
python3 -m DIMet.src.isotopolog_prop_stacked.py toy_metabologram/analysis001/config-2-001.yml
python3 -m DIMet.src.differential_analysis toy_metabologram/analysis001/config-2-001.yml
python3 -m DIMet.src.metabologram toy_metabologram/analysis001/config_metabologram.yml
```
Or create your own .sh file to launch it in a smarter way.
We provide also other very easy-to-understand minimal examples
without metabolograms, to guide you about the type of inputs and analyses that are possible, 
see [examples/readme_examples.md](examples/readme_examples.md).

## Option 3 : With Snakemake

This is the preferred way to use DIMet, because the reproducibility, 
and re-usability are very high with snakefiles, and our tool is more 
suited to this . 

See :
- [tests/Snakefile_common.smk](tests/Snakefile_common.smk)
- [tests/Snakefile_with_metabologram](tests/Snakefile_with_metabologram)

## Option 4 : With Galaxy

If you feel more comfortable in the web, try our newest Galaxy : 

-----------------------------------


## Measures' files requirements

- xlsx files are not admitted. Check if Tracegroomer may help you with your case. 

- The Isotopologues' names  must be written with a '**\_m+**' separating the name of the metabolite and the number of labeled carbons, example: 'Citrate\_m+0', 'Citrate\_m+1', and so on. Note: if you used our Tracegroomer,  this nomenclature is set automatically.

- The Isotopologue proportions must be comprised between 0 and 1 (0.1, 0.7, etc). The same is required for mean enrichment or fractional contributions.

### Input files

1. The measures' files : in the "data/" folder
2. The metadata file : in the "data/" folder
3. The .yml configuration file : in the "analysis--/" folder

Regarding the measure's files, they consist of 4 tab-delimited .csv files, one by type of measure:

- metabolites' abundances ("abundance")
- labelled (carbon) fractional contributions or mean enrichment ("meanE\_or\_fracContrib")
- isotopologues' proportions ("isotopologue\_prop")
- isotopologues' absolute values ("isotopologue\_abs")

You can also run DIMet if you have all except the last type of measure. 


Regarding the .yml file, we supply examples that you can use as template, such  [this config.yml template](examples/toy1/analysis001/config-1-001.yml). To double-check your modifications there exist online editors, such as https://yamlchecker.com/, just copy-paste the .yml file content, and easily edit!

  
        
 #### FAQ : 
 
 - I have my tracer metabolome measures in another very different format (a single xlsx, a single tsv, etc). Is there any fast method to generate the separated tables that DIMet requires?
 Yes, use https://github.com/johaGL/Tracegroomer
 
 
 - I have only 1 (or 2, or 3) of the 4 required measures', and I do have isotopologues' absolute values. Is there any fast method to generate the rest?
 Yes, use https://github.com/johaGL/Tracegroomer
 
 
 - I have only 3 measures' (because the software that I used for isotopologues' correction does not provide isotopologues' absolute values.) Can I still use DIMet ?
 Yes, you can use DIMet on the rest of the measures, check the [examples/toy2/](examples/toy2/)  and organize your job similarly.
 

 - Can I have only n=2 in each group ?
 Some of our "examples" contain only two biological replicates by timepoint and condition, for practical illustration purposes. Having only two replicates is discouraged (low statistical power). 


 - After the differential analysis, why do I obtain identical padj values ?
   The Benjamini-Hochberg formula uses the ranks of the p-values, when there are p-values that are very similar, it can happen that the computed _padj_ yields identical values, see : https://github.com/statsmodels/statsmodels/issues/4408   



## Contact us
Use the 'issues' of this github repo for any question that you could not answer by reading our guide.
Feel free to contact us so we can help you to make your data a suitable input for DIMet.



### References

[1] Bruntz, R. C., Lane, A. N., Higashi, R. M., & Fan, T. W. M. (2017). Exploring cancer metabolism using stable isotope-resolved metabolomics (SIRM). Journal of Biological Chemistry, 292(28), 11601-11609.

[2] Buescher, J. M. et al. A roadmap for interpreting (13)C metabolite labeling patterns from cells. Curr Opin Biotechnol 34, 189–201, https://doi.org/10.1016/j.copbio.2015.02.003 (2015).

[3] Guyon J, Fernandez‐Moncada I, Larrieu C, M, Bouchez C, L, Pagano Zottola AC, Galvis J,...& Daubon T, (2022). Lactate dehydrogenases promote glioblastoma growth and invasion via a metabolic symbiosis, EMBO Molecular Medicine, e15343

### Acknowledgments
To the VIB traning for the course in Tracer Metabolomics, 2023. 
To Thomas Daubon team and Joris Guyon for provided data, from we extracted toy examples.

