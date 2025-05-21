# behavior-epi-models
# behavior-epi-models

This file contains information for the code used to create all figures and carry out analyses in the publications [1]. MATLAB versions used were 2019a, 2023a, and 2024a. Descriptions of the code in each folder are below:

“Identifiability Code” contains ipynb and .nb files used in Theorems 3-6 (pg. 8-10, structural identifiability of the endogenous models) and Theorems B.1-B.4 (pg. 25 in Appendix B, structural identifiability of the exogenous models), which were run using the Julia package StructuralIdentifiability.jl [2].

“Figure 2 - Numerical Simulations” contains .m files for the numerical simulation graphics in Figure 2 (pg. 10), which were run using MATLAB.

“Figure 3 - Exog. time-vary sens” contains .m files used in the Figure 3 graphs (pg. 12) of time-varying sensitivity analysis for the exogenous systems, which were run using MATLAB. 

“Figure 4 - Endog. time-vary sens” contains .m files used in the Figure 4 graphs (p. 13) of time-varying sensitivity analysis for the endogenous systems, which were run using MATLAB. 

“Figure 5 - PRCC” contains .m files and .mat files used in the Figure 5 graphs (p. 15) of global sensitivity analysis for all four systems, which were run using MATLAB. Files should be run following these steps:

	To generate bar graphs:
	1	Exogenous models: run ‘Bars_SEIR_SEIRS.m’
	⁃	Needs input files ‘LHS_PRCC_SEIR_output.mat’ and ‘LHS_PRCC_SEIRS_output.mat’
	2	Endogenous models: run ‘Bars_SEIRb_SEIRSb.m’
	⁃	Needs input files ‘LHS_PRCC_SEIRb_output.mat’ and ‘LHS_PRCC_SEIRSb_output.mat’

	To generate output files:
	1	SEIR: run ‘LHS_PRCC_SEIR.m’
	⁃	Needs input file ‘LHS_raw.mat’
	2	SEIRb: run ‘LHS_PRCC_SEIRb.m’
	⁃	Needs input file ‘LHS_raw.mat’
	3	SEIRS: run ‘LHS_PRCC_SEIRS.m’
	⁃	Needs input file ‘LHS_raw.mat’
	4	SEIRSb: run ‘LHS_PRCC_SEIRSb.m’
	⁃	Needs input file ‘LHS_scaled_sample.mat’

“Figure C.1 - Endog change params” contains .m files used in the Figure C.1 graphs (p. 30 in Appendix C) showing how the the infectious population size of the SEIRSb system changed in response to changes in parameters, which were run using MATLAB.

“Figures 6,7,D.1,D.2,D.3 - Model Validation” contains .m files and .xlsx files (for data) used in Figure 6 (p. 17), Figure 7 (pg. 19), Table D.1 (pg. 32 in Appendix D), Figures D.1 and D.2 (pg. 33 in Appendix D), Table D.2 (p. 34 in Appendix D), and Figure D.3 (pg. 35 in Appendix D). The .m files were run using MATLAB.

Citation
If you find this repository valuable for your research, we kindly request that you acknowledge our paper by citing it.

[1] LeJeune L, Ghaffarzadegan N, Childs LM, Saucedo O. Mathematical analysis of simple behavioral epidemic models. Mathematical Biosciences. 2024 Sep 1;375:109250.

[2] Dong R, Goodbrake C, Harrington HA, Pogudin G. Differential elimination for dynamical models
via projections with applications to structural identifiability. SIAM Journal on Applied Algebra and
Geometry. 2023;7(1):194-235.
