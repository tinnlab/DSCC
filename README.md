# DSCC

This repo contains the source code and scripts for reproducing the results reported in our paper "DSCC: Disease subtyping using Spectral clustering and Community detection from Consensus networks".
Please follow the below instructions:

1. Requirements

You need an environment with both R and Python installed.
Also, each of the provided scripts might require some packages being installed before you can run them.
   
2. Setup
```bash
# Clone the repository
git clone https://github.com/tinnlab/DSCC.git
cd DSCC

# create a data folder
mkdir Data
cd /Data

# Download processed data
# Download data for DSCC
wget "https://seafile.tinnguyen-lab.com/f/9aa42f40d78247bc95db/?dl=1" -O DSCC_Main.zip
wget "https://seafile.tinnguyen-lab.com/f/2ad884778b00467ea690/?dl=1" -O DSCC_Relevant.zip

# Download data for comparison methods
wget "https://seafile.tinnguyen-lab.com/f/7c943f8328b3452db479/?dl=1" -O Others_Main.zip
wget "https://seafile.tinnguyen-lab.com/f/5f1d638d32a641d08e56/?dl=1" -O Others_Relevant.zip

# unzip the downloaded files
```

3. Run analysis for DSCC and comparisons method
```bash
# Please check the scripts for required packages before running them
# Run DSCC
Rscript Run_DSCC.R --no-save

# Run Comparison methods
Rscript Run_Other.R --no-save
```

4. Calculate the performance metrics
```bash
# Cox p-values and numbers of clusters
Rscript GetCoxPv.R --no-save

# C-Indices
Rscript GetData_SP.R --no-save
Rscript trainpredict_all_SP.R --no-save
python Eval_Methods.py
```
