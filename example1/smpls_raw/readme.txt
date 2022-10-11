Raw data 'MCF001089_6668_Cycloserine-metabolomics_results.xlsx'
 needs transformation because xlsx sheets have
samples in rows and metabolites in columns
and we need the opposite.

Steps:

1. titles were deleted, and a new sheet "samples_sheet" was added,
yielding : modifiedfile_Cycloserine.xlsx

2. The script formatter_be.py
using that file modifiedfile_Cycloserine.xlsx produced well
structured files, those  files were copied  in: example/data/ 
      
Note: 
- formatter_be.py is very data-dependent
    and is just an example of how to give correct formatting.
- formatter_be.py results have been moved to example/data/
- step 1  was done manually
    

