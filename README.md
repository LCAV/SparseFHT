SparseFHT Code repo
===================

This is the code for SparseFHT algorithm as presented in the following papers.

[Long] R. Scheibler, S.Haghighatshoar, and M. Vetterli,
[_A Fast Hadamard Transform for Signals with Sub-linear Sparsity in the Transform Domain_](http://arxiv.org/abs/1310.1803),
arXiv:1310.1803, 2013.

[Short] R. Scheibler, S. Haghighatshoar, and M. Vetterli,
[_A Fast Hadamard Transform for Signals with Sub-linear Sparsity_](http://infoscience.epfl.ch/record/189818),
Allerton conference 2013.


Version
-------

The code version is v0.2. Note that the code stability might not be optimal at
the moment. We are working on it!

Contact
-------

Robin Scheibler 
[email](mailto:robin.scheibler@epfl.ch)
[homepage](http://lcav.epfl.ch/Robin_Scheibler)

Please do not hesitate to contact me for help and support!
I would be happy to help you run the code.

Plateform
---------

The code has been tested on Mac OS X 10.7, 10.8, 10.9, and on Ubuntu linux.

Matlab mex wrappers were used for the code generating the figures in the paper. The core of the algorithm is implemented in C.

Compile
-------

    cd <code_dir> 
    make

To reproduce the figures from the paper, type in the following in a matlab shell:

    cd matlab/
    compile_c_code
    ErrorSim
    ErrorSim_plot
    TimingSim
    TimingSim_plot
    LessSparseSim
    LessSparseSim_plot

__Note__

1. The simulation is fairly time-consuming.
2. To speed-up the simulation, the parallel toolbox was used. If you do not have the parallel toolbox, replace the `parfor` instructions by `for` in all the Sim scripts.

License
-------

2013 (c) Robin Scheibler, Saeid Haghighatshoar, Martin, Vetterli.

The code is free to reuse for academic purposes. For any other purposes, contact the authors.
