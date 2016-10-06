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
             'data/20161005-164723_timing_sim.pickle',
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

    # check if a processed pickle file exists for these files
    pickle_file = os.path.splitext(data_files[0])[0] + '_proc_{}'.format(len(data_files)) + '.pickle'

    if os.path.isfile(pickle_file):
        # read the pickle file
        with open(pickle_file, 'rb') as f:
            df, parameters = pickle.load(f)
            f.close()

    else:
        # build the data table line by line
        print 'Building table...'
        columns = ['n', 'algorithm', 'time', 'k', 'C']
        table = []

        # This is the output from `figure_doa_experiment.py`
        for data_file in data_files:

            with open(data_file, 'rb') as f:
                data = pickle.load(f)
                f.close()

            # extra variables
            parameters = data['parameters']
            args = data['args']
            sim_out = data['out']

            for i,a in enumerate(args):
                k_lst = np.arange(1,a[0]-1)
                c_lst = np.zeros(k_lst.shape)
                print sim_out[i][0].shape
                print k_lst.shape
                for u in range(sim_out[i][0].shape[0]):  # loops
                    for v in range(sim_out[i][0].shape[1]):  # k value
                        table.append([a[0], 'SparseFHT', sim_out[i][0][u,v], k_lst[v], c_lst[v]])
                for t in sim_out[i][1]:  # loops
                    table.append([a[0], 'FHT', t, np.nan, np.nan])

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

    I = df['algorithm'] == 'FHT'
    df_fht = df[I][['n','time', 'algorithm']]
    df_fht = df_fht.pivot_table(values='time', columns='n', aggfunc=np.mean)

    df_sfht = df[df['algorithm'] == 'SparseFHT'].pivot_table(index='k', values='time', columns='n', aggfunc=np.mean)


    n = parameters['n']
    alpha_star = np.zeros(n.shape[0])

    # alpha_star is defined as the smallest
    for i,nn in enumerate(n):
        pos = np.where(df_sfht[nn] - df_fht[nn] > 0)[0][0]-1
        if pos < 0:
            alpha_star[i] = 0
        else:
            alpha_star[i] = df_sfht.index[pos] / nn


    # Plot random
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

    plt.plot(parameters['n'], alpha_star)

    # remove the x-grid
    ax = plt.gca()
    ax.xaxis.grid(False)

    plt.xlabel('$\\log_2N$')
    plt.ylabel('$\\alpha^\\star$')
    plt.xlim([parameters['n'][0], parameters['n'][-1]])
    plt.ylim([-0.05, 0.85])

    sns.despine(offset=10, trim=False, left=True, bottom=True)

    plt.tight_layout(pad=0.5)

    if not os.path.exists('./figures'):
        os.mkdir('./figures')

    plt.savefig('./figures/alpha_star.pdf', format='pdf')
    plt.savefig('./figures/alpha_star.png', format='png')

    # Second figure illustrate the runtime for a specific size
    plt.figure()

    nn = 15

    kk = np.arange(1, nn-1)
    plt.semilogy(kk, df_sfht[nn][:nn-2] / 1000)
    plt.semilogy(kk, np.ones(kk.shape)*df_fht[nn] / 1000)

    # remove the x-grid
    ax = plt.gca()
    ax.xaxis.grid(False)

    ax.text(4, 0.60*df_sfht[nn][4]/1000, 'SparseFHT')
    ax.text(4, 1.15*df_fht[nn]/1000, 'FHT')

    plt.xlabel('$\\alpha$')
    plt.ylabel('Runtime [s]')
    plt.xlim([kk[0], kk[-1]])
    plt.xticks([nn / 3, 2 * nn / 3], ['1/3', '2/3'])
    plt.ylim([0.95*np.min(df_sfht[nn]/1000), 1.05*np.max(df_sfht[nn]/1000)])

    sns.despine(offset=10, trim=False, left=True, bottom=True)

    plt.tight_layout(pad=0.5)

    plt.savefig('./figures/runtime.pdf', format='pdf')
    plt.savefig('./figures/runtime.png', format='png')
