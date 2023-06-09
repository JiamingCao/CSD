# CSD

Codes for generating results in paper "*A Model of Neurovascular Coupling and its Application to Cortical Spreading Depolarization*"

The whole model consists of three major parts: ion exchange model, neurovascular coupling model, and hemodynamic model. The first two parts are described by functions "genTuckwell2D" and "nvcoupling", respectively. The hemodynamic model is fairly straightforward without having to solve differential equations, so is simply a section in "RunThis.m"

List of files:

- genTuckwell2D(ratio)
  
  - "ratio" is the percentage power of neuronal reuptake [0, 1]

- nvcoupling(t, y, K_probe, Glu_probe, T, CPP, serca_ratio, BK_shift)
  
  - This is the system function of the neurovascular coupling model that must be supplied to one of the ODE solvers. Here we used the Matlab built-in ode15 function
  
  - t is the time vector, "y" is a place holder required by ODE solver functions, "K_probe""Glu_probe" are results from the previous section, "T" is the time vector of K_probe and Glu_probe, "CPP" is chosen to be 80mmHg for all the simulations, "serca_ratio" is the efficiency of the SERCA pump [0,1], and "BK_shift" is the negative shift of one of the BK channel parameters that governs its sensitivity

- vasodynamic(t, y, K, CPP)
  
  - This is a subsection of the nvcoupling model
  
  - K: perivascular potassium concentration, CPP: CPP value in mmHg

- sweepK
  
  - Simulates vasomotion in relationship with perivascular potassium concentration and CPP. Generates Fig. 8

- RunThis
  
  - The main code that runs all the models and generates all figures except for Fig. 8
  
  - Please pay attention to the comments as well as the manuscript that indicates which parameters are to be used for each figure
