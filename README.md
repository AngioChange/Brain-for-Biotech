Brain-for-Biotech (BfBio) is a graph-based tool built for gene prioritization. Based on a Personalized PageRank (PPR) [1,2] algorithm applied on an integrated network of different omics datasets and publicly available gene-gene/protein-protein interaction databases, it allows to infer biological functions in genes for which is was not known before based on a *guilt-by-association* principle. BfBio ranks all the genes in the biological network according to their connectivity with *seed genes*, i.e. genes already studied and characterized in a specific context.

[1] Brin S, Page L. The anatomy of a large-scale hypertextual Web search engine. *Computer Networks and ISDN Systems.* 1998 Apr;**30**(1):107-17. doi:10.1016/S0169-7552(98)00110-X.

[2] Page L, Brin S, Motwani R, Winograd T. The PageRank Citation Ranking: Bringing Order to the Web. *Stanford Digital Library Technologies Project.* 1998. doi:10.1.1.31.1768.

We describe BfBio in the following preprint: [BfBio: a graph-based tool for the prediction of Angiogenic Stalk Cell genes using a Personalized PageRank algorithm.](https://www.biorxiv.org/content/10.64898/2026.08.23.746493v3)

# Installation
In order to get BfBio, you will need to clone this repository:

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


