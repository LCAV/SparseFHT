from __future__ import division

import sys
import copy
import pickle
import numpy as np
import pandas as pd
import getopt
import os

import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == "__main__":

    argv = sys.argv[1:]
    data_files = [
             'data/20161005-005059_error_sim.npz',
             ]

    try:
        opts, args = getopt.getopt(argv, "hf:", ["file=",])
    except getopt.GetoptError:
        print('figure_doa_separation_plot.py -f <data_file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('figure_doa_separation_plot.py -f <data_file>')
            sys.exit()
        elif opt in ("-f", "--file"):
            data_files = arg.split(',')

    # check if a pickle file exists for these files
    pickle_file = os.path.splitext(data_files[0])[0] + '_{}'.format(len(data_files)) + '.pickle'

    if os.path.isfile(pickle_file):
        # read the pickle file
        with open(pickle_file, 'rb') as f:
            df, parameters = pickle.load(f)
            f.close()

    else:
        # build the data table line by line
        print 'Building table...'
        columns = ['K', 'B', 'C', 'algorithm', 'seed', 
                'MSE', 'supp size', 'bit error', 'success', 'unsat', 'loops']
        table = []

        # This is the output from `figure_doa_experiment.py`
        for data_file in data_files:
            data = np.load(data_file)

            # extra variables
            parameters = data['parameters'].tolist()
            args = data['args'].tolist()
            sim_out = data['out'].tolist()

            IL = parameters['inner_loops']
            for i,a in enumerate(args):
                for j in range(IL):
                    table.append(a + sim_out[i*IL + j])

        # create a pandas frame
        print 'Making PANDAS frame...'
        df = pd.DataFrame(table, columns=columns)

        # turns out all we need is the follow pivoted table
        #perf = pd.pivot_table(df, values='Error', index=['SNR'], columns=['Algorithm'], aggfunc=np.mean)

        with open(pickle_file, 'wb') as f:
            pickle.dump([df, parameters], f)
            f.close()

    sns.set(style='whitegrid')
    sns.plotting_context(context='poster', font_scale=2.)
    pal = sns.cubehelix_palette(8, start=0.5, rot=-.75)

    # Draw the figure
    print 'Plotting...'

    df_rand = df[df['algorithm'] == 'RANDOM']
    df_det = df[df['algorithm'] == 'DETERMINISTIC']

    # Plot random
    p_rand = pd.pivot_table(df_rand, values='success', index=['K'], columns=['C'], aggfunc=np.mean)
    p_rand = p_rand.reindex_axis(sorted(p_rand.columns, key=int), axis=1)
    p_rand = p_rand.reindex_axis(sorted(p_rand.index, key=int), axis=0)

    p_det = pd.pivot_table(df_det, values='success', index=['K'], columns=['C'], aggfunc=np.mean)
    p_det = p_det.reindex_axis(sorted(p_det.columns, key=int), axis=1)
    p_det = p_det.reindex_axis(sorted(p_det.index, key=int), axis=0)

    
    sns.set(style='white', context='paper', font_scale=1.2,
            rc={
                'figure.figsize':(7.,3.15), 
                'lines.linewidth':2.,
                'font.family': 'sans-serif',
                'font.sans-serif': [u'Helvetica'],
                'text.usetex': False,
                })
    #pal = sns.cubehelix_palette(6, start=0.5, rot=-0.75, dark=0.25, light=.75, reverse=True, hue=0.9)
    pal = sns.cubehelix_palette(6, start=0.5, rot=-0.5,dark=0.3, light=.75, reverse=True, hue=1.)
    cmap = sns.cubehelix_palette(6, start=0.5, rot=-0.5,dark=0.3, light=.75, reverse=True, hue=1., as_cmap=True)
    sns.set_palette(pal)
    #sns.set_palette('viridis')

    plt.figure()

    # the axis values
    k = np.log2(np.array([int(K) for K in p_rand.index]))
    C = np.array([int(c) for c in p_rand.columns])
    extent = [int(k[0]), int(k[-1]), int(C[0]), int(C[-1])]

    # Theoretical value of C
    kp = np.r_[k, k[-1]+1]
    n = parameters['n']
    cth = np.empty_like(kp)
    cth[:np.floor(n/3)] = n/kp[:np.floor(n/3)]
    cth[np.floor(n/3):np.floor(2*n/3)] = 3
    cth[np.floor(2*n/3):] = n / (n - kp[np.floor(2*n/3):])

    # deterministic
    plt.subplot(1,2,1)
    plt.imshow(p_rand.as_matrix().T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
    plt.plot(kp-0.5, cth, 'k')

    plt.xlim([k[0], k[-1]])
    plt.ylim([C[0], C[-1]])
    plt.xticks([int(parameters['n'] / 3), int(2 * parameters['n'] / 3)], ['1/3','2/3'])
    plt.xlabel('Sparsity $\\alpha$')
    plt.ylabel('Oversampling $C$')

    sns.despine(offset=10, trim=False, left=True, bottom=True)

    # random
    plt.subplot(1,2,2)
    plt.imshow(p_det.as_matrix().T, origin='lower', aspect='auto', extent=extent, cmap=cmap)
    plt.plot(kp-0.5, cth, 'k')

    plt.xlim([k[0], k[-1]])
    plt.ylim([C[0], C[-1]])
    plt.xticks([int(parameters['n'] / 3), int(2 * parameters['n'] / 3)], ['1/3','2/3'])
    plt.yticks([])
    plt.xlabel('Sparsity $\\alpha$')

    sns.despine(offset=10, trim=False, left=True, bottom=True)

    plt.tight_layout(pad=0.5)

    plt.show()
