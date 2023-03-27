# DIMet

[![License](https://img.shields.io/github/license/johaGL/dimet?label=license)](https://github.com/johaGL/dimet/blob/main/LICENSE)

[![CBiB Logo](imgs/cbib_logo.png)](https://www.cbib.u-bordeaux.fr/)

## DIMet: Differential Isotope-labeled targeted Metabolomics
----------------------------
**DIMet** is a bioinformatics pipeline for differential analysis of isotopic targeted labeling data.

Closely related to conventional metabolomics, stable isotope-resolved metabolomics (SIRM) uses an isotope labeled substrate to track specific pathways. From these data, it is possible to compute differences in isotope enrichment, changes in the labeling pattern, or differences in the contribution of nutrients to a metabolite pool, yielding knowledge of the metabolic state [1, 2]. Targeted metabolomics, when combined to transcriptomics, allows to better characterize perturbations within the pathways of interest.  

DIMet supports the analysis of full metabolite abundances and isotopologue contributions, and allows to perform it either in the differential comparison mode or as a time-series analysis. As input, the DIMet accepts three types of measures: a) isotopologues’ contributions, b) fractional contributions (also known as mean enrichment), c) full metabolites’ abundances. Specific functions process each of the three types of measures separately.

DIMet is intended for downstream analysis in corrected tracer data (corrected for the presence of natural isotopologues). Make sure you that the metabolomics platform provides you the output of the correction procedure before using our DIMet pipeline. 

![schema](imgs/schemaAlone4github.png)



##### Table of Contents  
[Requirements](#requirements)  

[Run the Analysis (fast guide)](#run-the-analysis-fast-guide)
- [Prepare](#prepare)  
- [Get PCA(s)](#get-pcas)
- [Run Differential analysis](#differential-analysis)
- [Get plots](#get-plots)
- [..](#get-grams)
- 
[Detailed guide](#detailed-guide)
- [Prepare in depth](#prepare-in-depth)


## Requirements

You need a UNIX system, with conda or miniconda3 installed, see [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Then clone DIMet repository in your `$HOME` directory, and set the virtual environment by running (from terminal):
```
cd $HOME
git clone git@github.com:cbib/DIMet.git
conda env create --file DIMet/dimet.yml
conda activate dimet 
```


# Run the Analysis (fast guide)


## Prepare


DIMet `prepare` module is the first step of the entire downstream analysis offered by our pipeline. It is required to be performed before any of the other modules. Its execution takes only few seconds! 

### Input files

Three types of Tracer data are accepted by our `prepare` module:
1. IsoCor results (.xlsx measurments file).
2. Results provided by the VIB Metabolomics Expertise Center (El-Maven results are shaped by VIB MEC team into a multi-sheet .xlsx file).  
3. A 'generic' .xlsx measurments file.

Your data is expected to correspond to one of the three types above. We provide you 'toy' examples of each one of them  [here](examples/readme_examples.md). Pick the example most suited to your data.

For example, let's take [toy3 example](examples/toy3/) , which contains:

 - 'data/' folder, with:
    * the tracer metabolomics data as one .xlsx file. 
    * the samples description that we call the "metadata", see below
    
 - 'analysis001/' folder, with:
     * your configuration file (extension .yml), see below.



Regarding the "metadata", we explain it in detail in the section [Metadata](#the-metadata).


Regarding the .yml file, we supply examples that you can use as template, such  [this one](examples/toy1/analysis001/config-1-001.yml). To double-check your modifications there exist online editors, such as https://yamlchecker.com/, just copy-paste and edit!


### Execute `prepare` 

The examples serve to demonstrate how fast this module can be. Take toy3 example, copy and paste the entire toy3 folder in your 'home/' folder, then from terminal:
```
python -m DIMet.src.prepare toy3/analysis01/config-3-01.yml
```


### Output files

The prepare output consist of 4 tables inside the folder: 'analysis001/results/prepared/tables/', one for each type of measurement, by compartment (for example, if you have 'cellular' and 'supernatant' compartments, there will be 8 tables). Each table contains the metabolites as rows, and the samples as columns. 

### Final words about `prepare`
Be patient, if this is the first time you perform downstream Tracer Metabolomics analysis with our DIMet pipeline, it can take some time to get used to prepare your input files. But once you have done it once, it will be really fast to analyze new datasets, do not give up. Check our section [Detailed guide](#detailed-guide). 

## Get PCA(s)

In construction

## Get plots

In construction

## Differential analysis

In construction

## Get grams

In construction

# Detailed guide

## Prepare in depth

To remind, your data is already the result of the procedure correcting the areas or intensities for the presence of isotopologues in nature. There are several software options (for example, IsoCor, El-Maven, etc) to perform that correction. After correction, you can use our DIMet for the downstream analysis. DIMet `prepare` module is the first step of this downstream analysis. 

In any of the three possible scenari (IsoCor, VIB or generic), you we encourage you to organize your files as shown in the examples that we provide. If will help you to have your results and your configuration files close, for reproductibility.





### Advanced prepare options

We provide advanced options for `prepare` module, check the help:
```
python -m DIMet.src.prepare --help
```
they appear as 'optional arguments' in the help menu.


You can:

- normalize by the amount of material (number of cells, tissue weight): setting the path to the file in your **.yml** configuration. The file must be like [this one](examples/toy2/data/nbcells-or-amountOfMaterial.csv), and the first column must contain the same names as in metadata 'former\_name'.
- normalize by an internal standard (present in your data) at choice: using the advanced option `--use_internal_standard`.

However we have some indications that can slightly differ for [users having VIB results as input](#users-having-vib-results), [users having IsoCor results](#users-having-isocor-results) or users having ['generic' type of data](#users-having-generic-data).

### Users having VIB results

As shown in the example 'toy2' [here](examples/toy2/), give the names of the sheets that are present in your excel file coherently. 
 
Our pipeline performs, by default:
- the subtraction of the means of the blanks across all metabolites' abundance for each sample.
- seting to NaN the values of abundance that are under the limit of detection (LOD).
- excluding metabolites whose abundance values across all samples are under LOD (excluded then from all tables by default).
- stomping fractions values to be comprised between 0 and 1 (some negative and some superior to 1 values can occur after )

You can modify all those options depending on your needs, they appear as 'optional arguments' in the help menu. 

 
### Users having IsoCor results

A typical IsoCor results table is described in: https://isocor.readthedocs.io/en/latest/tutorials.html#output-files
 
 Our pipeline transforms its columns into tables, so here the 'Isocor column : DIMet table' correspondence:
    - corrected_area : isotopologuesAbsolute  
    - isotopologue_fraction : isotopologuesProportions
    - mean\_enrichment :  meanEnrich\_or_fracContrib
    - (nothing)   : Abundance 

Abundance table is the sum (per metabolite) of IsotopologuesAbsolute, we also perform it automatically, as this column is not present in the input data.     
    
Please stick to the example [toy1](examples/toy1/) for the names of the tables in the **.yml** file for isocorOutput
    
        
Options regarding to detection limit (LOD) and blanks will not have any effect on the analysis. LOD is not provided in the data, and the same is true for blanks. 

All the other options do have effect: those related to internal standard, amount of material, and isotopologues.
 
 

### Users having generic data

We have created this option for those formats that are not the other two scenarios, so your data is expectd to be in the form of a .xlsx file with sheets similar as in VIB results:
- sheets corresponding to isotopologue Proportions (when available) and isotopologue Absolute values must have isotopologues as columns and samples as rows.
- sheets corresponding to abundance and mean enrichment  must have metabolites as columns and samples as rows.

As in example [toy3](examples/toy3) if you only have isotopologue Absolute values, but not the other tables: put them as a single named sheet in your .xlsx file, and we automatically generate all the other types of tables for you ! 



End of the remarks regarding the advanced options.


 
## The "Metadata"

Here the first lines of the required metadata table, which must be a .csv (comma delimited) file : 

| sample            | timepoint | condition | timenum | short_comp  |  former_name |
|-------------------|-----------|-----------|-------|------------|--------------- |
| Control\_cell\_T0-1 | T0        | Control   | 0     | cell       | MCF001089_TD01 |
| Control\_cell\_T0-2 | T0        | Control   | 0     | cell       | MCF001089_TD02 |
| Control\_cell\_T0-3 | T0        | Control   | 0     | cell       |  MCF001089_TD03|

You can create it with any spreadsheet program such as Excel or Google Sheets or LibreOfice Calc. At the moment of saving your file you specify that the delimiter must be a comma, see https://support.microsoft.com/en-us/office/save-a-workbook-to-text-format-txt-or-csv-3e9a9d6c-70da-4255-aa28-fcacf1f081e6. 

Column names in metadata must be exactly: 
 - former\_name
 - sample
 - timepoint
 - timenum
 - condition
 - short\_comp

 
 
The column 'former\_name' must have the names of the samples **as given in your data**. 
  
 
 The column 'sample' must have the names as you want them to be (or set identical to former\_name if you prefer). To set  names that are meaningful is a better choice, as we will take them for all the results.
 
 
 
 The column 'timenum' must contain only the numberic part of the timepoint, for example 2,0, 10, 100  (this means, without letters ("T", "t", "s", "h" etc) nor any other symbol). Make sure these time numbers are in the same units (but do not write  the units here!).
 
 

The column 'short\_comp' is an abbreviation, coined by you, for the compartments. This will be used for the results' files names: the longer the compartments names are, the longer the output files' names! Please pick short and clear abbreviations to fill this column.


## Contact us
Use the 'issues' of this github repo for any question that you could not answer by reading our guide.
Feel free to contact us so we can help you to make your data a suitable input for DIMet.


 
### References

[1] Bruntz, R. C., Lane, A. N., Higashi, R. M., & Fan, T. W. M. (2017). Exploring cancer metabolism using stable isotope-resolved metabolomics (SIRM). Journal of Biological Chemistry, 292(28), 11601-11609.

[2] Buescher, J. M. et al. A roadmap for interpreting (13)C metabolite labeling patterns from cells. Curr Opin Biotechnol 34, 189–201, https://doi.org/10.1016/j.copbio.2015.02.003 (2015).

[3] Guyon J, Fernandez‐Moncada I, Larrieu C, M, Bouchez C, L, Pagano Zottola AC, Galvis J,...& Daubon T, (2022). Lactate dehydrogenases promote glioblastoma growth and invasion via a metabolic symbiosis, EMBO Molecular Medicine, e15343





## things to do
- go to settings-> action -> enable repo modification

