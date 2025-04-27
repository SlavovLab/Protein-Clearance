# Protein-Clearance
Analysis from Leduc and Slavov, *bioRxiv*, doi: [10.1101/2025.02.10.637566](https://doi.org/10.1101/2025.02.10.637566): **Impact of protein degradation and cell growth on mammalian proteomes**

## Data
Protein clearance and synthesis rates were measured from murine tissue samples using *in vivo* metabolic pulse with lysine labeled with stable heavy isotopes. Samples were then analyzed by (i) timsTOF Ultra 2, which quantified over 12,000 total protein concentrations and over 10,000 degradation rates across all tissues and (ii) Exploris 480, which quantified 8,740 concentrations and degradation rates:

![](Protein-Degradation-MouseTissues-Datasets.png)

- Raw LC-MS and searched data can be found on MassIVE, [MSV000097050](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=732135590e454112b315f6e610891080).
- Searched and processed data can be found on [Zenodo](https://zenodo.org/records/14827610).

----- 
  

## Analysis 
In this work, we demonstrate that protein degradation plays a significant role in determining protein concentrations across the proteome in slowly growing cells. Reduction in the influence of degradation in growing cells and variation in rates of degradation across samples contribute to changes in protein concentration across different conditions. 

To reproduce these finding, the following two R scripts can be run:
- Tissue_metabolic_labeling.R for plotting figures based on tissue data aqcuired by Leduc *et al.*, 2025
- External_data_and_simulation.R for plotting figures based on previously published data


Both scripts call functions from the Functions_for_analysis.R script. 
The starting point should be making sure all the packages that are loaded in the begining of the functions script are installed.

These scripts read in data files that are store on [Zenodo](https://zenodo.org/records/14827610) so there should be no need to download any additional data.

