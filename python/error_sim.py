from __future__ import division
import numpy as np

n = 22

# Global simulation parameter
params = {
        'n': n,
        'b': np.arange(1, n-1),
        'C': np.arange(1, 13),
        'N': 2**n,
        'algo_names': ['RANDOM', 'DETERMINISTIC'],
        'max_iter': 20,
        'sigma2': 100,
        'outer_loops': 100,
        'inner_loops': 10,
        'seed': 54321,
        }

# initialize RNG
np.random.seed(params['seed'])

# Prepare the arguments
args = []
for bval in params['b']:
    for Cval in params['C']:
        for algo in params['algo_names']:
            for i in range(params['outer_loops']):
                seed = np.random.randint(4294967295, dtype=np.uint32)
                args.append([2**bval, 2**bval, Cval, algo, seed])

# The local loop
def parallel_loop(args):

    K = args[0]
    B = args[1]
    C = args[2]
    algo_name = args[3]
    seed = args[4]

    import numpy as np
    import scipy.linalg as la
    import time

    import pysparsefht
    from utils import random_k_sparse

    try:
        import mkl as mkl_service
        # for such parallel processing, it is better 
        # to deactivate multithreading in mkl
        mkl_service.set_num_threads(1)
    except ImportError:
        pass


    # for such parallel processing, it is better 
    # to deactivate multithreading in mkl
    mkl_service.set_num_threads(1)

    if algo_name == 'RANDOM':
        algo = pysparsefht.ALGO_RANDOM
    elif algo_name == 'DETERMINISTIC':
        algo = pysparsefht.ALGO_OPTIMIZED
    else:
        ValueError('No such algorithm.')

    # initialize rng
    np.random.seed(seed)

    # a list for return values
    ret = []

    # Run the inner loops
    for i in range(params['inner_loops']):

        # generate a seed for the C RNG
        C_seed = np.random.randint(4294967295, dtype=np.uint32)

        # create sparse vector
        x_hat, y_hat, supp_hat = random_k_sparse(params['N'], K, params['sigma2'])

        # compute WHT
        x = pysparsefht.fht(x_hat)

        # Now apply the SparseFHT
        y_hat2, supp_hat2, unsat, loops = pysparsefht.sparse_fht(x, int(K), int(B), int(C),
                                    max_iter=params['max_iter'],
                                    algo=algo,
                                    req_loops=True, req_unsat=True,
                                    seed=C_seed)

        supp_r = supp_hat2[supp_hat2 > 0]
        y_r = y_hat2[supp_hat2 > 0]

        # reconstructed vector
        x_hat2 = np.zeros(params['N'])
        x_hat2[supp_r] = y_r

        mse = la.norm(x_hat - x_hat2)**2
        supp_size = supp_r.shape[0]
        bit_error = len(set(supp_r).symmetric_difference(set(supp_hat)))
        success = (unsat[0] == 0)

        ret.append([mse, supp_size, bit_error, success, unsat, loops])

    return ret


if __name__ == '__main__':

    '''
    This script will run a monte-carlo simulation to determine the probability
    of failure of the SparseFHT for different parameters of the algorithm
    '''

    import os, sys, getopt
    import numpy as np
    import scipy.linalg as la
    import time

    import pysparsefht


    # default values
    serial_flag = False
    test_flag = False
    data_filename = None

    # parse arguments
    cmd_name = sys.argv[0]
    argv = sys.argv[1:]

    def print_help(cmd):
        print('%s [-t -s] -f <filename>' % cmd)
        print('  -s, --serial: Use serial computing')
        print('  -t, --test: Test mode (run 1 loop)')
        print('  -f <filename>, --file=<filename>: name of output file')

    try:
        opts, arguments = getopt.getopt(argv, "hf:ts", ["file=", "test","plot"])
    except getopt.GetoptError:
        print_help(cmd_name)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print_help(cmd_name)
            sys.exit()
        elif opt in ("-f", "--file"):
            data_filename = arg
        elif opt in ("-t", "--test"):
            test_flag = True
        elif opt in ("-s", "--serial"):
            serial_flag = True

    #------------------#
    # LOCAL PARAMETERS #
    #------------------#

    # There is the option to only run one loop for test
    if test_flag:
        print 'Running one test loop only.'
        args = []
        for bval in params['b']:
            for Cval in params['C']:
                for algo in params['algo_names']:
                    seed = np.random.randint(4294967295, dtype=np.uint32)
                    args.append([2**bval, 2**bval, Cval, algo, seed])

    # Main processing loop
    if serial_flag:
        print 'Running everything in a serial loop.'
        # Serial processing
        out = []
        for ag in args:
            out.append(parallel_loop(ag))

    else:
        import ipyparallel as ip

        print 'Using ipyparallel processing.'

        # Start the parallel processing
        c = ip.Client()
        NC = len(c.ids)
        print NC,'workers on the job'

        # Push the global config to the workers
        c[:].push(dict(params=params))

        # evaluate the runtime
        then = time.time()
        out1 = c[:].map_sync(parallel_loop, args[:NC])
        now = time.time()
        one_loop = now - then
        print 'Total estimated processing time:', len(args)*one_loop / len(c[:])

        # dispatch to workers
        out = out1 + c[:].map_sync(parallel_loop, args[NC:])

    # unwrap the outer and inner loops
    out_unrolled = [inner for outer in out for inner in outer]

    # Save the result to a file
    if data_filename is None:
        date = time.strftime("%Y%m%d-%H%M%S")
        data_filename = 'data/{}_error_sim.npz'.format(date)
    
    np.savez(data_filename, args=args, parameters=params, out=out_unrolled) 

    print 'Saved data to file: ' + data_filename
