# LSTM-GLM_LSTM-NN
This repository contains the Python and R code used to conduct the simulation studies presented in
"Dynamic Prediction of Residual Life with Longitudinal Covariates Using Long Short-Term Memory Networks" 
by Grace Rhodes, Marie Davidian, and Wenbin Lu.


We provide 8 programs, which should be executed in numerical order according to their file names.
A description of the 8 programs can be found below.

1) 01_generate_data.R (R program)
   creates the simulated data sets of patient covariates and survival times.

2) 02_create_contextVecs.py (Python Program)
   creates the window-specific context vectors for each simulation.

3) 03_create_linReg_mixedEff.R (R Program)
   creates the linear regression vectors and the mixed effects vectors 
   for each simulation.

4) 04_create_fpca.R (R Program)
   creates the FPCA vectors for each simulation.

5) 05_trainTest_ipcw.R (R Program)
   divides the patients in each simulation into training and testing data sets
   and computes the inverse probability of censoring weight of each patient.

6) 06_lstmNN_pred.py (Python Program)
   fits the LSTM-NN for each simulation and outputs the predicted MRLs.

7) 07_loss_cIndex.R (R Program)
   computes the training and testing losses and c-indexes
   of the LSTM-GLM, the LSTM-NN, and the six dynamic transformed MRL models 
   (baseline, last-value carried forward, average, linear regression, mixed effects, FPCA)
   for each simulation.

8) 08_create_boxplots.R (R Program)
    creates the boxplots of the training and testing losses and c-indexes. 

To facilitate running the aforementioned programs, we recommend creating the following directory structure:

Create a parent directory, "Simulations." 

Create the following sub-directories in "Simulations": data, tuner, context_vectors, linReg_mixedEff, fpca, ipcw, lstmNN, results.

In each sub-directory, create the sub-sub-directories: aft, cox.
