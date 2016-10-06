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
             'data/20161005-115611_less_sparse_sim.npz',
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
        columns = ['alpha', 'seed', 
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

        with open(pickle_file, 'wb') as f:
            pickle.dump([df, parameters], f)
            f.close()

    sns.set(style='whitegrid')
    sns.plotting_context(context='poster', font_scale=2.)
    pal = sns.cubehelix_palette(8, start=0.5, rot=-.75)

    # Draw the figure
    print 'Plotting...'

    # Plot random
    pf = df.groupby(['alpha']).mean()
    pf = pf.reindex_axis(sorted(pf.index, key=float), axis=0)
    
    sns.set(style='whitegrid', context='paper', font_scale=1.2,
            rc={
                'figure.figsize':(4.,3.), 
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
    beta = parameters

    plt.plot(parameters['beta'], pf['success'].as_matrix())

    ax = plt.gca()

    # remove the x-grid
    ax.xaxis.grid(False)


    plt.xlabel('$\\beta$')
    plt.ylabel('Success probability')

    plt.xticks([1., 2., 3., 4.])
    plt.yticks([0., 0.25, 0.5, 0.75, 1.])
    plt.xlim([0.99, 4.01])
    plt.ylim([-0.01, 1.01])

    sns.despine(offset=10, trim=False, left=True, bottom=True)

    plt.tight_layout(pad=0.5)

    if not os.path.exists('./figures'):
        os.mkdir('./figures')

    plt.savefig('./figures/less_sparse.pdf', format='pdf')
    plt.savefig('./figures/less_sparse.png', format='png')
