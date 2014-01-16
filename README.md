SparseFHT Code repo
===================

This is the code for SparseFHT algorithm as presented in the following papers.

[Long] R. Scheibler, S.Haghighatshoar, and M. Vetterli,
[_A Fast Hadamard Transform for Signals with Sub-linear Sparsity in the Transform Domain_](http://arxiv.org/abs/1310.1803),
arXiv:1310.1803, 2013.

[Short] R. Scheibler, S. Haghighatshoar, and M. Vetterli,
[_A Fast Hadamard Transform for Signals with Sub-linear Sparsity_](http://infoscience.epfl.ch/record/189818),
Allerton Conference on Communication, Control and Computing, 2013.

Abstract
--------

A new iterative low complexity algorithm has been presented for computing the
Walsh-Hadamard transform (WHT) of an N dimensional signal with a K-sparse WHT,
where N is a power of two and K = O(N^α), scales sub-linearly in N for some 0 <
α < 1. Assuming a random support model for the non- zero transform domain
components, the algorithm reconstructs the WHT of the signal with a sample
complexity O(K log\_2(N/K)), a computational complexity O(K log\_2(K) log\_2(N/K))
and with a very high probability asymptotically tending to 1.

The approach is based on the subsampling (aliasing) property of the WHT, where
by a carefully designed subsampling of the time domain signal, one can induce a
suitable aliasing pattern in the transform domain. By treating the aliasing
patterns as parity- check constraints and borrowing ideas from erasure
correcting sparse-graph codes, the recovery of the non-zero spectral values has
been formulated as a belief propagation (BP) algorithm (peeling decoding) over
a sparse-graph code for the binary erasure channel (BEC). Tools from coding
theory are used to analyze the asymptotic performance of the algorithm in the
“very sparse” (α ∈ (0,1]) and the “less sparse” (α ∈ (1,1)) regime.

Contact
-------

Robin Scheibler 
[email](mailto:robin[dot]scheibler[at]epfl[dot]ch)
[homepage](http://lcav.epfl.ch/Robin_Scheibler)

Please do not hesitate to contact me for help and support!
I would be happy to help you run the code.

Plateform
---------

The code has been tested on Mac OS X 10.7, 10.8, 10.9, and on Ubuntu linux.

Matlab mex wrappers were used for the code generating the figures in the paper. The core of the algorithm is implemented in C.

Run the code in Matlab
----------------------

### Reproduce the figures from the paper

To reproduce the figures from the paper, type in the following in a matlab shell:

    cd <path_to_SparseFHT>/matlab/
    make_mex_files
    make_figures

__Note__

1. The simulation is fairly time-consuming.
2. To speed-up the simulation, the parallel toolbox was used. If you do not have the parallel toolbox, replace the `parfor` instructions by `for` in all the Sim scripts.

### Mex files info

* **SparseFHT** The fast sparse Hadamard transform algorithm.

         [Y, S, U, I] = SparseFHT(X, K, B, C, L, T)

         Wrapper for the Sparse Fast Hadamard Transform
         
         Input arguments:
         X: input vector  (size n)
         K: the sparsity (and size of y)
         B: number of buckets
         C: oversampling factor
         L: maximum number of iterations of decoder
         T: Type of algorithm to use ('Random' / 'Deterministic' / 'Optimized')
         
         Output arguments: 
         Y: output vector (size K)
         S: support vector (size K)
         U: the number of unsatisfied checks (optional)
         I: the number of loops run (optional)


* **FastHadamard** A straighforward implementation of the conventional fast
  Hadamard transform. No scaling factor is applied, i.e. applying the algorithm
  twice will result in the input vector weighted by its length.

        [Y] = FastHadamard(X)

        Input arguments:
        X: input vector (size has to be a power of two)

        Output arguments:
        Y: output vector, the Hadamard transform of X.

* **HadamardBenchmark** Calls a C routine performing a timing comparison of SparseFHT and FastHadamard.

        [Tfht Tsfht] = HadamardBenchmark(N, K, B, C, L, R, SEED)

        Input arguments:
        N: transform size to investigate (scalar, power of two)
        K: sparsity parameter, number of non-zero tranform domain coefficients (vector, power of two)
        B: number of buckets to use in SparseFHT (vector, same size as K, power of two)
        C: oversampling factor
        L: maximum number of iterations of decoder
        R: a length 4 vector containing the following parameters
            1. Number of repetitions of one measurement
            2. Number of warm-up run 
            3. Number of iterations for one measurement
            4. Maximum magnitude of non-zero components in the sparse signal
        SEED: A seed for the C random number generator

        Output arguments:

        Tfht: The runtime measurement of FastHadamard
        Tsfht: An array containing the runtime measurement of SparseFHT for every value of K

Reuse the C code
----------------

The makefile in `./C` folder will compile all the code as well as a number
of example/test files. This can be used as a basis to reuse the C code directly.

    cd ../C
    make all

The libraries included are:

* SparseFHT,
* Fast Hadamard tranform,
* Fast linear algebra in GF2 (boolean matrices and vectors algebra),
* Test files.

License
-------

2013-2014 (c) Robin Scheibler, Saeid Haghighatshoar, Martin, Vetterli.

The code is free to reuse for academic purposes. However, please acknowledge its use with the following citation.

    @InProceedings{EPFL-CONF-189818,
       author      = {Scheibler, Robin and Haghighatshoar, Saeid and Vetterli, Martin},
       booktitle   = {51{s}t {A}nnual {A}llerton {C}onference on
                     {C}ommunication, {C}ontrol, and {C}omputing},
       location    = {Allerton Retreat Center, Monticello, Illinois},
       title       = {A {F}ast {H}adamard {T}ransform for {S}ignals with
                     {S}ub-linear {S}parsity},
       year        = 2013,
       ee          = {http://arxiv.org/abs/1310.1803}
    }

For any other purposes, please contact the authors.
