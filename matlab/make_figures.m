% This file will run all the numerical experiments and reproduce the
% figures in the paper:
% R. Scheibler, S. Haghighatshoar, M. Vetterli,
% A Fast Hadamard Transform for Signals with Sub-linear Sparsity in the Transform Domain

% The numerical experiment scripts are split into two files.
% The first script runs the experiment and saves all the results in a mat file.
% The second script plots the data and format the figures.
% This allows to run first the experiments on a server or such, and later
% load the mat file and plot the data.

% Fig. 5 - The density evolution curve
DensityEvolution

% Fig. 7 - Error probability with deterministic hashing
Repetitions = 1000;
Type = 'Deterministic';
ErrorSim;
ErrorSim_plot;

% Fig. 8 - Error probability with random hashing
Repetitions = 1000;
Type = 'Random';
ErrorSim;
ErrorSim_plot;

% Fig. 9 - Threshold behavior in less-sparse regime (deterministic hashing)
Repetitions = 1000;
Type = 'Deterministic';
LessSparseSim;
LessSparseSim_plot;

% Fig. 10 and 11 - Runtime measurements
Repetitions = 1000;
TimingSim;
TimingSim_plot;
