Brain-for-Biotech (BfBio) is a graph-based tool built for gene prioritization. Based on a Personalized PageRank (PPR) [1,2] algorithm applied on an integrated network of different omics datasets and publicly available gene-gene/protein-protein interaction databases, it allows to infer biological functions in genes for which is was not known before based on a *guilt-by-association* principle. BfBio ranks all the genes in the biological network according to their connectivity with *seed genes*, i.e. genes already studied and characterized in a specific context.

[1] Brin S, Page L. The anatomy of a large-scale hypertextual Web search engine. *Computer Networks and ISDN Systems.* 1998 Apr;**30**(1):107-17. doi:10.1016/S0169-7552(98)00110-X.

[2] Page L, Brin S, Motwani R, Winograd T. The PageRank Citation Ranking: Bringing Order to the Web. *Stanford Digital Library Technologies Project.* 1998. doi:10.1.1.31.1768.

We describe BfBio in the following preprint: [BfBio: a graph-based tool for the prediction of Angiogenic Stalk Cell genes using a Personalized PageRank algorithm.](https://www.biorxiv.org/content/10.64898/2026.08.23.746493v3)

# Installation
In order to install BfBio, you will need to clone this repository:

```
git clone https://github.com/AngioChange/Brain-for-Biotech
cd Brain-for-Biotech
```
The easiest way to prepare the environment with all required packages is to use [conda](docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). You can create the proper environment with:

```
conda env create -n BfBio -f BfBio.yaml
```
Remember to activate your environment before running the scripts and notebooks!
```
conda activate BfBio
```

# Main functionalities of BfBio
## Performing gene prioritization
BfBio's purpose is to predict the functional roles of under-annotated genes by extracting patterns from complex biological networks. Even though the networks described in the preprint have been built to identify genes important for vascular endothelial cells (EC), but BfBio's framework can be applied to any other specific biological process as long as relevant co-expression datasets, protein-protein interaction (PPI) databases and training genes are available. 

## Fetching relevant Pubmed publications through text mining
Although most human protein coding genes have functional annotations in databases, such as [GeneCards](www.genecards.org), many remain poorly characterized. Considering there is a daunting and constant need for new drugs, these poorly characterized genes represent a potential goldmine for knowledge gain and drug discovery. To investigate whether the predicted targets are underexplored or already well defined in the context of interest, we have written a script, based on the existing [Biopython](www.biopython.org) library, to look for relevant PubMed publications associated with the predicted genes. 

## Notebooks
The notebooks available from this GitHub describe several frameworks, such as BfBio's main workflow, the preprocessing of the annotation-based databases considered in this work, or several benchmarking procedures. Feel free to adapt them to your needs. 

