# Library
import os
from tqdm import tqdm 

# Listing all the correlation threshold values that will be generated
corr_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

for corr in tqdm(corr_values, desc = "Processing correlation thresholds ..."):
	cmd = f"conda activate BfBio ; $VSC_SCRATCH/projects/BfBio_reloaded/scripts/Rscript Network_Analysis4Omics.R {corr} inflammation graph"
	os.system(cmd)