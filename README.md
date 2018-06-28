# bayNorm: relevant code for producing figures in the paper
code for producing figures in bayNorm


#Purpose of this repository
The main purpose of this repository is to provide the analysis procedure used in the paper.


# Source code of bayNorm
Source code of bayNorm can be found [here](https://github.com/WT215/bayNorm)


# Real datasets used in this paper
This paper involves the following 8 studies:

1. Klein study (https://www.cell.com/cell/abstract/S0092-8674%2815%2900500-0)
2. Grün study (https://www.nature.com/articles/nmeth.2930)
3. Torre study (https://www.cell.com/cell-systems/abstract/S2405-4712(18)30051-6)
4. Bacher study (https://www.nature.com/articles/nmeth.4263)
5. Islam study (https://www.ncbi.nlm.nih.gov/pubmed/21543516)
6. Soumillon study (https://www.biorxiv.org/content/early/2014/03/05/003236)
7. Tung study (https://www.nature.com/articles/srep39921)
8. Patel study (http://science.sciencemag.org/content/344/6190/1396)

# Simulated datasets used in this paper:
There are 4 simulated datasets with DE genes. Each one of them consists of 2 two groups of cells, and 100 cells in each group. 2000 out of 10000 genes were simulated to be DE genes in the first group and half of the 2000 genes were upregulated `\Simulations\SIM_DE`.

1. SIM DE I: mean capture efficiency $<\beta>=10\%$ for two groups.
2. SIM DE II: mean capture efficiency $<\beta>=5\% \text{ and } 10\%$ for two groups respectively.
3. SIM DE III: mean capture efficiency $<\beta>=10\% \text{ and } 5\%$ for two groups respectively.
4. SIM DE IV: mean capture efficiency $<\beta>=5\% \text{ and } 5\%$ for two groups respectively.


There are another 2 simulated datasets without DE genes `\Simulations\SIM_noDE`. Mean capture efficiency $<\beta>=10\% \text{ and } 5\%$ for two groups respectively. These two simulations were inspired by Bacher study. The purpose is to study the ability of normalization method in terms of correcting different sequencing depths.

1. SIM Bacher I: Parameters were estimated from Klein study.
2. SIM Bacher II: Parameters were estimated from H1_P24 cells from Bacher study.

# Datasets and the corresponding figures

## Real datasets
1. Klein study: Fig1 (b)-(e), Fig3 (a)-(b); Fig S2, S8a-b.
2. Grün study: Fig2 a,c,e and g; Fig S11a-b, S12-S13
3. Torre study: Fig2 b,d,f and h; Fig S6, S8e-f, S10a, S11c, S14.
4. Bacher study: Fig S7, S9 a-d, S10e, S16, S19a, S23a.
5. Islam study: Fig 3c; Fig S9e-f, S23b.
6. Soumillon study: Fig3d; Fig S21.
7. Tung study: Fig4, FigS3-S5, S8c-d, S10b-d, S25-26
8. Patel study: FigS10f

## Simulated datasets with DE genes

1. SIM DE I: FigS15a,e,i, S20c-d, S22a, S24a, S27-29
2. SIM DE II: FigS15b,f,j, S20c-d, S22b, S24b, S27-29
3. SIM DE III: FigS15c,g,k, S20c-d, S22c, S24c, S27-29
4. SIM DE IV: FigS15d,h,l, S20c-d, S22d, S24d, S27-29

## Simulated datasets without DE genes
1. SIM Bacher I: S17, S19b, S20a-b
2. SIM Bacher II: S18, S19c


## Some notes before running the code
1. **You cannot directly run all the code at the same time. The paths in each R file need to be modified accordingly.** 
2. The normalization and DE detection could take a long time, which depends on the size of raw data. Hence make sure running the code step by step so as to avoid bugs. 
3. Useful functions are stored in the file `\Functions`, some of them need to be loaded in advance. 
4. The noramlization method `DCA` is developed using python. The Jupyter Notebooks for running DCA are stored in the file `\DCA`. Make sure running DCA normalization and corresponding DE detection, and them feed the DCA normalized data into the other R files.  
5. Some R files need several `.RData` files as input and will also output `.RData` files used in other cases. Hence make sure the first step is completed so as to produce necessary `.RData` files to begin with. 







# The first step

## Preparing for the real datasets

1. Klein study: firstly, run `\RealData\Klein_study\Klein_bayNorm.R`, output `Klein_bayNorm.RData`.

2.Grün study: run `LOAD_Grun_smFISH.R` (output `smFISH_norm_load.RData`), `LOAD_Grun_2i.R` (output `Grun_2014_RAW.RData`) and `LOAD_Grun_serum.R` (output `Grun_2014_RAW_serum.RData`). Then run `Grun_2i_norms.R` (output `Grun_2i_norms.RData`) and `Grun_serum_norms.R` (output `Grun_serum_norms.RData`) for normalizing data. Note that the other method DCA needs to be run separately.

3. Torre study: run `Load_Torre.R` (output `Load_Torre.RData`). Then run `Torre_many_normalizations.R` (out put `Torre_many_normalizations.RData`) for normalizing data. 

4. Bacher study: run `LOAD_Bacher.R` (output `RAW_INITIATE.RData`) to load H1 and H9 datasets. Then run `H1_many_normalizations.R` (output `"H1_many_normalizations.RData"`) and `H9_many_normalizations.R` (output `"H9_many_normalizations.RData"`) respectively.


5. Islam study: run `Load_Islam.R` (output `Load_Islam.RData`). Then run `Islam_many_normalizations.R` (output `Islam_many_normalizations.RData`).
6. Soumillon study: run `LOAD_Soumillon.R` (output `Soumillon_2014.RData`). Then run `Soumillon_norms.R` (output `Soumillon_analysis.RData`).
7. Tung study: run `Load_Tung.R` (output `Load_Tung.RData`). Then run `Tung_many_normalizations.R` (output `Tung_norms.RData`).
8. Patel study: run `Load_Patel.R` (output `Patel2014_bay_out.RData`)


## Notes before running simulations
Firstly, we need to estimate parameters from the real data. Relevant codes are stored in `\bayNorm_papercode\Figure1`.

1.  For Klein dataset, if you have completed the step 1 as shown above, then `Klein_bayNorm.RData` stored the parameters you need. `Klein_bayNorm.RData` is needed in SIM DE I-IV and  SIM Bacher I.
2. For Bacher dataset (H1_P24), run a section named `REAL DATA 6: Bacher study (H1_P24 cells)` in the file `Simulations_realdata.R`, which output `H1p24_bay_sim_allgene.RData` used in SIM Bacher II.






## Preparing for the simulated datasets (with DE genes)
The codes are stored in: `\Simulations\SIM_DE`
1. SIM DE I: run `DE_sim_01_01.R` (output `SIM_1.RData` and `GG_SIM_1.RData`).
2. SIM DE II: run `SIM_005_01.R` (output `SIM_005_01.RData` and `GG_SIM_005_01.RData`).
3. SIM DE III: run `SIM_01_005.R`(output `SIM_01_005.RData` and `GG_SIM_01_005.RData`).
4. SIM DE IV: run `SIM_005_005.r`(output `SIM_005_005.RData` and `GG_SIM_005_005.RData`).

## Preparing for the simulated datasets (without DE genes)
The codes are stored in: `\Simulations\SIM_noDE`
1. SIM Bacher I: run `SIM_noDE_01_005.R` (output `SIM_noDE_01_005.RData`)
2. SIM Bacher II: run `SIM_noDE_01_005_H1.R` (output `SIM_noDE_01_005_H1.RData`)










