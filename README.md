# bayNorm: relevant code for producing figures in the paper
code for producing figures in bayNorm


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

# The first step: preparing for the raw data












