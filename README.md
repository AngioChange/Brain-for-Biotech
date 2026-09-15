Brain-for-Biotech (BfBio) is a graph-based tool built for gene prioritization. Based on a Personalized PageRank (PPR) [1,2] algorithm applied on an integrated network of different omics datasets and publicly available gene-gene/protein-protein interaction databases, it allows to infer biological functions in genes for which is was not known before based on a *guilt-by-association* principle. BfBio ranks all the genes in the biological network according to their connectivity with *seed genes*, i.e. genes already studied and characterized in a specific context.

[1] Brin S, Page L. The anatomy of a large-scale hypertextual Web search engine. *Computer Networks and ISDN Systems.* 1998 Apr;**30**(1):107-17. doi:10.1016/S0169-7552(98)00110-X.
[2] Page L, Brin S, Motwani R, Winograd T. The PageRank Citation Ranking: Bringing Order to the Web. *Stanford Digital Library Technologies Project.* 1998. doi:10.1.1.31.1768.

# Installation
In order to get BfBio, you will need to clone this repository:

```
git clone https://github.com/AngioChange/Brain-for-Biotech
cd Brain-for-Biotech
```
The easiest way to prepare the environment with all required packages is to use conda. You can create the proper environment with:

```
conda env create -n BfBio -f BfBio.yaml
```
Remember to activate your environment before running the scripts and notebooks!
```
source active BfBio
```


