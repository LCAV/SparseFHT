"""
Pysparsefht module provides a nice wrapper around
the C functions computing conventional and sparse Hadamard transforms.

Example
-------

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

"""

import numpy as np
import sparsefht_wrapper

def fht(x):
    '''
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
    '''

    if x.ndim == 1:
        n = x.shape[0]
    elif x.ndim == 2:
        n = x.shape[1]
    else:
        raise ValueError("The vector x should be 1D or 2D.")

    if 2**int(np.log2(n)) != n:
        raise ValueError("The transform length should be a power of two.")

    if x.dtype != np.double:
        x = np.array(x, dtype=np.double)

    return sparsefht_wrapper.fht(x, np.empty_like(x))

def sparse_fht(x, K, B=None, C=3, max_iter=20, algo=sparsefht_wrapper.ALGO_OPTIMIZED, req_loops=False, req_unsat=False, seed=0):
    '''
    Wrapper for the Sparse Fast Hadamard Transform

    Parameters
    ----------
    x: ndarray (1D or 2D)
        Input vector, the size should be a power of two.
    K: int
        The sparsity expected
    B: int, optional
        The number of buckets, by default the closest power of two larger than K is used
    C: int, optional
        The oversampling factor, (default: 3)
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
        The output vector of locations of non-zero coefficients (size K), if less coefficients than K
        were found, the remaining values are padded with -1.
    unsat: int or ndarray
        The number of unsatisfied checks (if req_unsat == True)
    loops: int or ndarray
        The number of loops run by decoder (if req_loops == True)
    '''

    if B is None:
        B = int(2**np.ceil(np.log2(K)))
        print B
        print C
        print algo

    if x.ndim == 1:
        n = x.shape[0]
        num_transforms = 1
        out_shape = (K)
    elif x.ndim == 2:
        n = x.shape[1]
        num_transforms = x.shape[0]
        out_shape = (num_transforms, K)
    else:
        raise ValueError("The vector x should be 1D or 2D.")

    if 2**int(np.log2(n)) != n:
        raise ValueError("The transform length should be a power of two.")

    if 2**int(np.log2(B)) != B:
        raise ValueError("The number of buckets should be a power of two.")

    # check type of x
    if x.dtype != np.double:
        x = np.array(x, dtype=np.double)

    # create output vectors
    y = np.zeros(out_shape, dtype=np.double)
    support = np.zeros(out_shape, dtype=np.int32)

    if req_loops:
        loops = np.zeros(num_transforms, dtype=np.int32)
    else:
        loops = np.zeros(0, dtype=np.int32)

    if req_unsat:
        unsat = np.zeros(num_transforms, dtype=np.int32)
    else:
        unsat = np.zeros(0, dtype=np.int32)

    # Call the C code
    sparsefht_wrapper.sparse_fht(x, y, support, unsat, loops, 
                                int(K), int(B), int(C), 
                                int(max_iter), int(algo), int(seed))

    if not (req_unsat or req_loops):
        return y, support
    elif req_unsat and not req_loops:
        return y, support, unsat
    elif not req_unsat and req_loops:
        return y, support, loops
    else:
        return y, support, unsat, loops

def benchmark(sparsities, buckets, oversampling, transform_size, 
        loops=10, warm=0, body=1, max_mag=500, sfht_max_iter=20, seed=0):
    '''
    Run a benchmark to compare the runtime of sparse_fht vs conventional fht

    Parameters
    ----------
    sparsities: ndarray
        A vector that contains a number of values for the sparsity
    buckets: ndarray
        A vector that contains the number of buckets to use for a specific sparsity
    oversampling: int or ndarray
        A scalar or vector of values for the oversampling factor of the sparse_fht
    transform_size: int
        The transform size considered
    loops: int, optional
        The number of loops to run (default: 10)
    warm: int, optional
        The number of warm-up run of the algorithm to get the machine in steady state (default: 0)
    body: int, optional 
        The number of repetitions to run for one timing measurement (default: 1)
    max_mag: int, optional
        The maximum magnitude of the sparse coefficients (default: 500)
    sfht_max_iter: int, optional
        The maximum number of iterations of the iterative decoder
    seed: int, optional
        A seed to pass to the random number generator (default: 0)

    Returns
    -------
    Tsfht: ndarray
        An array containing the timing value for each loop for the sparse_fht (loops x sparsities.shape[0])
    Tfht: ndarray
        An array containing the timing value for each loop for the fht (length loops)
    '''

    oversampling = np.array(oversampling, dtype=np.int32)
    if oversampling.ndim == 0:
        oversampling = np.array([oversampling])

    if sparsities.ndim != 1 or buckets.ndim != 1:
        ValueError("Benchmark: sparsities and buckets should be vectors")

    if sparsities.shape[0] != buckets.shape[0]:
        ValueError("Benchmark: the sparsity and bucket vectors should be the same length.")
    if sparsities.shape[0] != oversampling.shape[0] and oversampling.shape[0] != 1:
        ValueError("Benchmark: oversampling should be a scalar, or a vector the size of sparsity")

    if 2**int(np.log2(transform_size)) != transform_size:
        raise ValueError("The transform size should be a power of two.")

    if sparsities.dtype != np.int32:
        sparsities = np.array(sparsities, dtype=np.int32)
    if buckets.dtype != np.int32:
        buckets = np.array(buckets, dtype=np.int32)

    # The benchmark parameter array
    params = np.array([loops, warm, body, max_mag], dtype=np.int32)

    # The output arrays
    Tfht = np.zeros(loops, dtype=np.double)
    Tsfht = np.zeros((loops, sparsities.shape[0]), dtype=np.double)

    # Run the benchmark
    sparsefht_wrapper.benchmark(sparsities, buckets, oversampling, params, Tsfht, Tfht, 
            int(transform_size), int(sfht_max_iter), int(seed))

    return Tsfht, Tfht
    
