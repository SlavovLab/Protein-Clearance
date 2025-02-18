# Protein-Clearance
Analysis from Leduc and Slavov, *bioRxiv*, doi: [10.1101/2025.02.10.637566](https://doi.org/10.1101/2025.02.10.637566): **Protein degradation and growth dependent dilution substantially shape mammalian proteomes**

- Raw LC-MS and searched data can be found on MassIVE, [MSV000097050](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=732135590e454112b315f6e610891080).
- Searched and processed data can be found on [Zenodo](https://zenodo.org/records/14827610).

In this work, we demonstrate that protein degradation plays a significant role in determining protein concentrations across the proteome in slowly growing cells. Reduction in the influence of degradation in growing cells and variation in rates of degradation across samples contribute to changes in protein concentration across different conditions. 

To reproduce these finding, the following two R scripts can be run:
- External_data_and_simulation.R for plots from figure 1 and 2
- Tissue_metabolic_labeling.R for plots from figure 3 and 4

Both scripts call functions from the Functions_for_analysis.R script. 
The starting point should be making sure all the packages that are loaded in the begining of the functions script are installed.

These scripts read in data files that are store on [Zenodo](https://zenodo.org/records/14827610) so there should be no need to download any additional data.

