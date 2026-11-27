# AF3PoolSorter
Sort a proteome into pools for whole-proteome interaction screen with one query protein/protein complex.

## Description

A sorting algorithm that distributes all proteins from a FASTA file into pools containing a protein of interest, to do proteome-wide AlphaFold interaction screens.  

The output is a folder with json files, each containing up to 30 jobs for the AlphaFold3 server, with one pool each. The query can be a multi-, heteromeric protein complex.

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
The arguments should be:  
	1. The proteome fasta file (file ending in .txt, .fasta or .faa)  
	2. The ID of the protein of interest. Quotation marks are not necessary but they may save you trouble if there are special characters  
	3. The number of copies of the bait protein to include in the pool. If missing defaults to 1.  
	
### Heteromultimer query
You can specify any list of proteins to include in the query complex by passing them as an argument to the script. You can then specify the copy number of each protein in the complex by passing them too, in the order of the proteins specified. For example:

```
Rscript AF3PoolSorter.R NA1000proteome.faa YP_002518697.1  YP_002515378.1 2 3 

```

will generate pools containing 2 copies of YP\_002518697.1 and 3 copies of YP\_002515378.1. The order in which they are written will be the order of the chains in the output.


## Output
The code generates a folder containing all the json files. Since the protein IDs are not included in the pool files, a csv table is generated to map the IDs to the proteins in each pool.  

If there are proteins in the proteome so big that they wouldn't fit in the 5000 aa limit with the query protein, they are excluded from the search and a list of them is created.
## Considerations

AlphaFold3 limits the job size to 5000 aa. The algorithm takes a few steps to ensure that each pool is within the limit, but this is not mathematically guaranteed. If any pool were to exceed the size limit, a warning is printed on the terminal, as well as the identity of the pool. If this happens, you can manually remove a sequence from said pool in the json file (or contact me, maybe I can fix it).

# ipTM exporter
Pull ipTM numbers from a zip folder of AlphaFold results from jobs from AF3PoolSorter.

## Use
On AlphaFold server, select all jobs belonging to one screen and download in a single zip folder. The name of the jobs should have the form "Pool [number]".  

Move the zip folder and the script 'ipTM_exporter.R' into your working directory, as well as the csv map generated with AF3PoolSorter. You then pass those files as arguments to the script  as follows:

```
Rscript ipTM_exporter.R folds_2026_11_21_02_49.zip example_map.csv
```

## Output

A folder is generated with the name of the query protein. In the folder there's the following files
### All_ipTMs.csv
A very big table containing the scores for all pairwise combinations of the proteins in every pool. The values are the average of the top 5 models produced by AlphaFold for each pool.
**Chain 1 & 2**: Protein identifier within pool. A corresponds to the query protein, potentially B, C... too if the submission was a multimer  
**Corrected ipTM**: pairwise ipTM score, corrected by substracting the expected ipTM value, as described in Todor et al. 2026 (more on that below)  
**Protein 1 & 2**: Name of the proteins in the pair  
**Length 1 & 2**: Length of the proteins in the pair  
**Pair length**: Total length of the proteins in the pair, used to calculate the expected ipTM  
**Pool**: Pool identifier  
**Pair ipTM**: raw ipTM as calculated by AF3, averaged across the top 5 models for this pool  
**Expected ipTM**: Expected ipTM due to the size of complex. Exp ipTM = -0.036255571 + 0.004470512 * sqrt(Pair length)  
**ipTM standard deviation**: standard deviation of ipTM across the top 5 models for this pool  
**Pair pae min**: min pae per pair, averaged across the top 5 models for this pool  
**ranking score**: average ranking score from top 5 models of this pool  

### YourProtein.csv
The same table but only including pairs with the query protein as the first chain, sorted by corrected ipTM.

### ipTM correction plot
A plot of all pair length vs raw ipTM score, overlayed with the equation described by Todor et al. 2026, to get a sense of how accurate the correction is in this proteome.
