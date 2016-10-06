SparseFHT Code repo
===================

This is the code for SparseFHT algorithm as presented in the following papers.
The code is hosted on [github](https://github.com/LCAV/SparseFHT).

[Long] R. Scheibler, S.Haghighatshoar, and M. Vetterli,
[_A Fast Hadamard Transform for Signals with Sub-linear Sparsity in the Transform Domain_](http://infoscience.epfl.ch/record/204991),
IEEE Trans. Inf. Theory, vol. 61, 2015.

[Short] R. Scheibler, S. Haghighatshoar, and M. Vetterli,
[_A Fast Hadamard Transform for Signals with Sub-linear Sparsity_](http://infoscience.epfl.ch/record/189818),
Allerton Conference on Communication, Control and Computing, 2013.


Abstract
--------

In this paper, we design a new iterative low-complexity algorithm for
computing the Walsh-Hadamard transform (WHT) of an N dimensional signal with
a K-sparse WHT.  We suppose that N is a power of two and K = O(N^α), scales
sub-linearly in N for some α ∈ (0,1). Assuming a random support model for the
nonzero transform-domain components, our algorithm reconstructs the WHT of the
signal with a sample complexity O(K log\_(N/K)) and a computational complexity
O(K log\_2(K)log\_2(N/K)). Moreover, the algorithm succeeds with a high
probability approaching 1 for large dimension N.  

Our approach is mainly based on the subsampling (aliasing) property of the WHT,
where by a carefully designed subsampling of the time-domain signal, a suitable
aliasing pattern is induced in the transform domain. We treat the resulting
aliasing patterns as parity-check constraints and represent them by a bipartite
graph.  We analyze the properties of the resulting bipartite graphs and borrow
ideas from codes defined over sparse bipartite graphs to formulate the recovery
of the nonzero spectral values as a peeling decoding algorithm for a specific
sparse-graph code transmitted over a binary erasure channel (BEC). This enables
us to use tools from coding theory (belief-propagation analysis) to
characterize the asymptotic performance of our algorithm in the very sparse (α
∈ (0,1/3]) and the less sparse (α ∈ (1/3,1)) regime. Comprehensive simulation
results are provided to assess the empirical performance of the proposed
algorithm.

Contact
-------

Robin Scheibler 
([email](mailto:robin[dot]scheibler[at]epfl[dot]ch))
([homepage](http://lcav.epfl.ch/Robin_Scheibler))

Please do not hesitate to contact me for help and support!
I would be happy to help you run the code.

Plateform
---------

The code has been tested on Mac OS X 10.7, 10.8, 10.9, and on Ubuntu linux.

Python and Matlab mex wrappers were used for the code generating the figures in
the paper. The core of the algorithm is implemented in C.

The code needs a C99 compatible compiler, which seems not to be the case for
the windows version of Matlab currently. There seems to be some
[workaround](https://stackoverflow.com/questions/3737356/mex-problem-how-to-support-c99matlab),
one possibility being to rename all the `.c` files into `.cpp` and use a C++
compiler. All of this is untested as of now.

Run the code in Python
----------------------

A python wrapper for the C code was written and is available in `python/pysparsefht`.
The wrapper must be compiled as follows

    cd python/pysparsefht
    python setup.py build_ext --inplace

### Dependencies

The code relies on `numpy`, `scipy`, `matplotlib`, `seaborn`, `pandas`, and
`ipyparallel` for the parallel simulation. It is also possible to run all the
code serially, but this would take a long time.

### Reproduce the figures from the paper

A script is provided to reproduce all the figures in the paper.
Assuming the python module has been built as explained above, do the following.

    cd python/
    make_all_figures.sh

And that's it pretty much. The simulation scripts are suffixed with `_sim` and produce
a data file stored in `data` that can be reused later to generate the figure using
a second script suffixed with `_plot`. The figures will be stored in the `figures` folder.

The `make_all_figures.sh` scripts has a few options for testing and number of cores used.

    ./make_all_figures.sh [OPTS]
    Options:
      -t    Runs a single loop only for test purpose
      -s    Runs all the code in a simple for loop. No parallelism
      -n x  Runs the loops in parallel using x workers. This option is ignored if -s is used

The number of workers is in general set to the number of threads available minus one. That is
twice the number of cores, minus one.

Running `./make_all_figures.sh -t -n 7` on an Intel Core i7 2.8 GHz took 5 minutes.

### Pysparsefht module

A user-friendly python module was written around the C code and can be used to run the FHT and SparseFHT.

    import numpy as np
    import pysparsefht

    # Create a sparse Hadamard domain vector
    y = np.zeros(512)
    y[[13,72,121, 384]] = [12., -123., 91.5, -37.]

    # We can generate the dense time domain by
    # applying conventional FHT
    x = pysparsefht.fht(y)

    # Then, we use the sparse FHT to get back the
    # original sparse vector
    # y_val, y_loc contain the magnitude and locations, respectively, of the sparse signal
    y_val, y_loc = pysparsefht.sparse_fht(xs, 4)

    # The output of the transform is not normalize,
    # so we need to divide by the square root of the length
    y_hat = np.zeros(512)
    y_hat[y_loc] = y_val / np.sqrt(512)

    # says True
    np.allclose(y_hat, y)

The docstring for `sparse_fht`:

    Signature: pysparsefht.sparse_fht(x, K, B, C, max_iter=20, algo=3, req_loops=False, req_unsat=False, seed=0)
    Docstring:
    Wrapper for the Sparse Fast Hadamard Transform

    Parameters
    ----------
    x: ndarray (1D or 2D)
        Input vector, the size should be a power of two.
        K: int
          The sparsity expected
        B: int
          The number of buckets
        C: int
          The oversampling factor
        max_iter: int, optional
          The maximum number of iterations of decoder
        algo: int, optional
          The variant of the algorithm to use:
            * ALGO_RANDOM : Uses random hash functions
            * ALGO_DETERMINISTIC: Uses deterministic hash functions
            * ALGO_OPTIMIZED: Uses deterministic hash functions with optimized implementation (default)
        req_loops: bool, optional
          Requests to return the number of loops executed by the decoded
        req_unsat: bool, optional
          Requests to return the number of unsatisfied check nodes when the algorithm terminates
        seed: unsigned int, optional
          A seed for the random number generator to pass to the C code

    Returns
    -------
    y: ndarray (1D or 2D)
      The output vector of magnitudes (size K)
    support: ndarray (1D or 2D)
      The output vector of locations of non-zero coefficients (size K)
    unsat: int or ndarray
      The number of unsatisfied checks (if req_unsat == True)
    loops: int or ndarray
      The number of loops run by decoder (if req_loops == True)

The docstring for `fht`:

    Signature: pysparsefht.fht(x)
    Docstring:
    Fast Hadamard Transform

    This is a wrapper that calls a C implementation of the fast hadamard transform

    If x is a 2D array, the transform operates on the rows.

    The length of the transform should be a power of two.

    Parameters
    ----------
    x: ndarray (1D or 2D)
        The data to transform

        Returns
        -------
        A new ndarray containing the Hadamard transform of x




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

2013-2015 (c) Robin Scheibler, Saeid Haghighatshoar, Martin, Vetterli.

This work is licensed under the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy
of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

The code is free to reuse for non-commercial and academic purposes. However,
please acknowledge its use with the following citation.

    @article{EPFL-JOUR-204991,
       author      = {Scheibler, Robin and Haghighatshoar, Saeid and Vetterli, Martin},
       title       = {A {F}ast {H}adamard {T}ransform for {S}ignals with
                     {S}ub-linear {S}parsity in the {T}ransform {D}omain},
       journal     = {IEEE Trans. Inf. Theory}
       volume      = 61,
       year        = 2015,
       ee          = {http://infoscience.epfl.ch/record/204991}
    }

For any other purposes, please contact the authors.
