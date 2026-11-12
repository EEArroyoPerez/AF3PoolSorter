# AF3PoolSorter
Sort a proteome into pools for whole-proteome interaction screen with one query protein

## Description

A sorting algorithm that distributes all proteins from a FASTA file into pools containing a protein of interest, to do proteome-wide AlphaFold interaction screens.  

The output is a folder with json files, each containing up to 30 jobs for the AlphaFold3 server, with one pool each. The number of copies of the protein of interest can be controled  at the input, in case of multimeric proteins.

## Requirements
### Software
R version 4.6.1, with `dplyr`, `tidyr`, `stringr`, `seqinr` & `jsonlite` installed.
If not found, the packages will be installed on the first run of the script.
### Inputs 
FASTA file of the proteome.  
Identifier of the protein of interest, as used in the FASTA file.  

## Workflow

Copy the script 'AF3PoolSorter.R' into your working directory, as well as the proteome file.  
In the working directory, run the script as follows

```
Rscript AF3PoolSorter.R NA1000proteome.faa 'YP_002518697.1' 2
```
The arguments should be, in order:  
	1. The proteome fasta file  
	2. The ID of the protein of interest, with quotation marks  
	3. The number of copies of the protein to include in the pool. If missing defaults to 1.  
	
## Output
The code generates a folder containing all the json files. Since the protein IDs are not included in the pool files, a csv table is generated to map the IDs to the proteins in each pool.

## Considerations

AlphaFold3 limits the job size to 5000 aa. The algorithm takes a few steps to ensure that each pool is within the limit, but this is not mathematically guaranteed. If any pool were to exceed the size limit, a warning is printed on the terminal, as well as the identity of the pool. If this happens, you can manually remove a sequence from said pool in the json file (or contact me, maybe I can fix it).




