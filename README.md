# Hawkes process tutorials
The folder contains two tutorials on how to implement and use Bayesian Hawkes process models with the R-pacakge inlabru. Both of them use the Epidemic-Type-Aftershocks-Sequence (ETAS) model as example.

The folder how_to_build_Hawkes contains a tutorial describing how to implement Hawkes process models with the R-pacakge inlabru. It provides details on the functions to be provided by the user and how to pass them to inlabru. Also, the tutorial describes possible problems and difficulties and how to overcome them. It is based on the implementation of a temporal Epidemic-Type-Aftershock-Sequence model.

The folder how_to_use_Hawkes contains a tutorial describing how to use the functions contructed in the previous one. Specifically, I have wrapped the functions contained in how_to_build_Hawkes to facilitate their usage for non-expert users. I have provided one-line functions to read a seismic catalog, fit an ETAS model, retrieve the posterior distribution of the parameters, retrieve the posterior distribution of the number of events, obtain catalog-based forecasts for specified time intervals. All the funtions takes their input from a txt file which should be edited by the user. We provided two input file examples for the Amatrice and L'Aquila seismic sequences.

The folder code_for_paper contains all the code and data needed to reproduce the analysis presented in the paper at https://arxiv.org/abs/2206.13360
