# ABTCK

Augmented Bayesian treed co-kriging


## MIT License

Copyright (C) 2019 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

This file is part of the Github repository, ABTCK, Augmented Bayesian treed co-kriging

## Description 

MATLAB scripts implementing our proposed Bayesian procedure,'Augmented Bayesian treed co-kriging'. This procedure can be used for the statistical analysis of computer models when the available simulations are in multiple fidelity, they are not necessarily based on hierarchically nested experimental designs, and they may present local features, non-stationarities, or discontinuities.

This is a MATLAB code has been used to produce the numerical results and the statistical analysis reported in the paper  

* "Bayesian analysis of multifidelity computer models with local features and non-nested experimental designs: Application to the WRF model", by Alexandros Konomi and Georgios karagiannis; submited  in Journal of American Statistical Association for publication. 

## Requirements 

MATLAB compiler (R2017a or later)

## Files 

* cov_functions: Covariance functions neccessary for the GP 

* function_fold: Functions used in the simulation study

* Multi_GP_operations: MCMC operations for the co-kriging Gaussian process hyperparamter updates

* Multi_prediction: Functions used for prediction

* Multi_tree_nested: MCMC treed update of the model when the design is nested

* Multi_tree_NON_nested: MCMC treed update of the model when the design is nested

* Random_var: Various algorithms for generating data from distributions necessary in the BTC

* help_tree_operations: Operations that help to define the treed form

* SIMULATION_jasa: Contain the example 1 of the paper. It also contain different variations of that example. The user can choose to change the functiuon. 
This is not a fixed example -- you generate every time new data and check the performance of the algorithm. The performance may depend also on the simulated data. 

* Simulation_jasa2: Contain the example 2 of the paper. Here the user can also change the function and the setting
