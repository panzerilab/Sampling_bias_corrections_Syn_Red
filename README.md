# README

This directory contains the code and data used for the analyses in our paper. Below is a brief overview of the folder structure and key scripts:

## Folder Structure

- `Simulations_discrete.m`: MATLAB script for generating all discrete simulation datasets (except Figures 1, S12).
- `Simulations_Gaussian.py`: Python script for generating all Gaussian simulation datasets (except Figures 1, S12, and S13).
- `PlotFigures.m`: MATLAB script for plotting all simulation figures (and generating results for Figures 1 and S13).
- `tools.py`: collection of Python functions used by the `Simulations_Gaussian.py` script.
- `Figure_S12/`: Scripts to generate data and plots for Figure S12.
- `Real_data_analyses/`: Raw data and scripts for real data preprocessing and analysis.
- `Results/`: Stores all simulated results.
- `Figures/`: Stores all generated figures.
- `Functions/`: Helper functions, including implementations of our bias-correction algorithms. Python functions used in `Simulations_Gaussian.py` are included in `tools.py`.
  
## Setup Instructions


- To set up the Python environment, use the provided `environment.yml` file. If you're using Conda, run the following command to create the environment automatically:
  ```bash
  conda env create -f environment.yml
  ```
- Install the [MINT toolbox](https://github.com/panzerilab/MINT) for discrete PID computations in MATLAB.
- Install the [gPID toolbox](https://github.com/praveenv253/gpid) for Gaussian PID computations in Python.
- Install the [CaImAn toolbox](https://github.com/flatironinstitute/CaImAn) for the processing of the spike data for figure 3A
## Generating Spike Train Data

### For Figure 6A (Mouse Auditory Cortex)

1. Download data from: https://drum.lib.umd.edu/items/30d43732-7149-4726-a860-0ae3d210b2ae  
2. Add CaImAn to the MATLAB path:
   ```matlab
   addpath(genpath('CaImAn-MATLAB-master/'))
    ```
3. Generate the spike train data:
    ```
    data_for_information_analysis;
    generate_spike_train_data;
    ```
This will create the `spike_train_alltrials` variable saved in `DFF_spktrain_alltrials.mat`.

### Spike Trains Data for Figure 6B (Mouse Hippocampus)
1. Download the data from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012934, listed as 'S1_Data'
2. In the terminal, decompress and extract the archive:
   ```bash
   gunzip journal.pcbi.1012934.s002.gz
   tar -xvf journal.pcbi.1012934.s002
   ```
3.  Move the file `CA1_Data.mat` into `Real_data_analyses/Data`

### fMRI Data for Figure 6D
1. Download Q3 release from: https://www.humanconnectome.org/
2. Save it in `Real_data_analyses/Data`
