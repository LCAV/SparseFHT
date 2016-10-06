import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

if __name__ == "__main__":

    C = 3
    beta = np.array([3, 2, 1])
    pj = np.linspace(0, 1, 101)

    curves = (1 - np.exp(-beta[:,None] * pj[None,:]))**(C - 1)

    sns.set(style='white', context='paper', font_scale=1.2,
            rc={
                'figure.figsize':(3.15,3.15), 
                'lines.linewidth':2.,
                'font.family': 'sans-serif',
                'font.sans-serif': [u'Helvetica'],
                'text.usetex': False,
                })
            #pal = sns.cubehelix_palette(6, start=0.5, rot=-0.75, dark=0.25, light=.75, reverse=True, hue=0.9)
    pal = sns.cubehelix_palette(6, start=0.5, rot=-0.5,dark=0.3, light=.75, reverse=True, hue=1.)
    sns.set_palette(pal)
    #sns.set_palette('viridis')

    plt.figure()

    plt.plot(pj, pj, 'k--')

    for curve in curves:
        plt.plot(pj, curve)

    plt.xlabel('$p_j$')
    plt.ylabel('$p_{j+1}$')
    plt.xticks(np.linspace(0, 1, 5))
    plt.yticks(np.linspace(0, 1, 5))
    plt.axis('equal')

    ax = plt.gca()

    offset = [0.05, -0.12, -0.12]
    for b, curve, off in zip(beta, curves, offset):
        ax.text(0.65, curve[75] + off, '$\\beta = ' + str(b) + '$')

    sns.despine(offset=10, trim=False, left=False, bottom=False)

    plt.tight_layout(pad=0.5)

    if not os.path.exists('./figures'):
        os.mkdir('./figures')

    plt.savefig('./figures/density_evolution.pdf', format='pdf')
    plt.savefig('./figures/density_evolution.png', format='png')

