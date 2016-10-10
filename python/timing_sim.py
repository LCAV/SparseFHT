from __future__ import division
import numpy as np
import pickle

#---------------------------#
# Begin Configuration Zone  #
#---------------------------#

n = np.arange(5,24)   # Transform size 2**n

# Global simulation parameter
params = {
        'n': n,
        'algo_name': 'DETERMINISTIC',
        'max_iter': 10,
        'warm': 2,
        'body': 10,
        'max_mag': 500,
        'outer_loops': 100,
        'inner_loops': 10,
        'seed': 67584,
        }

# initialize RNG
np.random.seed(params['seed'])

# Generate the arguments
args = []
for nn in params['n']:
    for epoch in range(params['outer_loops']):
        seed = np.random.randint(4294967295, dtype=np.uint32)
        args.append([nn, seed])

#---------------------------#
# End of Configuration Zone #
#---------------------------#

# The local loop
def parallel_loop(args):

    import numpy as np
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

    n = args[0]

    b = np.arange(1, n-1)
    K = 2**b
    B = 2**b

    # compute value of C
    C = np.empty_like(b)
    C[:np.floor(n/3)] = n/b[:np.floor(n/3)]
    C[np.floor(n/3):np.floor(2*n/3)] = 3
    C[np.floor(2*n/3):] = n / (n - b[np.floor(2*n/3):])

    algo_name = params['algo_name']
    seed = args[1]

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

    # generate a seed for the C RNG
    C_seed = np.random.randint(4294967295, dtype=np.uint32)

    # create sparse vector
    Tsfht, Tfht = pysparsefht.benchmark(K, B, C, 2**n, 
            loops=params['inner_loops'], warm=params['warm'], body=params['body'], max_mag=params['max_mag'],
            sfht_max_iter=params['max_iter'], seed=C_seed)

    return [Tsfht, Tfht, b, C]


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
        params['inner_loops'] = 1
        params['warm'] = 0
        params['body'] = 1
        args = []
        for nn in params['n']:
            seed = np.random.randint(4294967295, dtype=np.uint32)
            args.append([nn, seed])

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

    # Save the result to a file
    if data_filename is None:
        date = time.strftime("%Y%m%d-%H%M%S")
        data_filename = 'data/{}_timing_sim.pickle'.format(date)
    
    # dump to pickle file
    with open(data_filename, 'wb') as f:
        pickle.dump(dict(args=args, parameters=params, out=out), f)
        f.close()

    print 'Saved data to file: ' + data_filename
