import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from utopya import DataManager, UniverseGroup
from utopya.plotting import is_plot_func, UniversePlotCreator
from .distribution_info import get_dist_info

# -----------------------------------------------------------------------------

#@is_plot_func(creator_type=UniversePlotCreator)
def media(dm: DataManager, *, out_path: str,
            uni: UniverseGroup):
    """
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (UniverseGroup): The selected universe data
    """
    #get data
    uni_cfg = uni['cfg']

    if uni_cfg['OpDyn']['media_status']=='on':
        num_steps = uni_cfg['num_steps']
        write_every = uni_cfg['write_every']
        num_users = uni_cfg['OpDyn']['nw_u']['num_vertices']
        avg_user_degree = uni_cfg['OpDyn']['nw_u']['mean_degree']
        num_media = uni_cfg['OpDyn']['nw_m']['num_vertices']
        data_m = uni['data/OpDyn/nw_media/opinion_m']
        user_numbers = uni['data/OpDyn/nw_media/user_count']
        time = np.linspace(0, num_steps/write_every, int(num_steps/write_every)+1)

        #set up format
        fig = plt.figure(figsize=(10,20), constrained_layout=True)
        widths = [6]
        heights = [0.5, 2, 3, 2]
        gs = fig.add_gridspec(ncols=len(widths), nrows=len(heights), width_ratios=widths, height_ratios=heights)

        #colour palette
        colors = {'media opinion':'dodgerblue'}
        fonts = {'title': 60, 'subtitle': 30, 'text':25, 'plot title': 20, 'label title': 18}

        #axis 0 is just for text
        axs0 = fig.add_subplot(gs[0, :])
        axs0.axis('off')
        axs0.text(0, 1.6, "Media parameter evolution", fontweight='bold', fontsize=30, verticalalignment='top', horizontalalignment = 'left')
        axs0.text(0, 1.0, "Evolution of media opinion and user numbers", fontsize=20)
        axs0.text(0, 0.4, "Numerical steps: "+str(num_steps), fontweight='bold', fontsize=13)
        axs0.text(0, 0.2, "Number of users: "+str(num_users), fontweight='bold', fontsize=13)
        axs0.text(0., 0., "Number of media: "+str(num_media), fontweight='bold', fontsize=13)

        #plot initial opinion state  .................................................................................................................................................
        axs1 = fig.add_subplot(gs[1, 0])
        axs1.set_title("Initial opinion distribution: "+get_dist_info(uni_cfg, 'media', 'opinion'), fontsize=11,  horizontalalignment = 'left', x=0)
        axs1.scatter(data_m[0,:], user_numbers[0, :], s=50, color=colors['media opinion'])
        axs1.xaxis.grid(lw=0.5)
        axs1.yaxis.grid(lw=0.5)
        axs1.set_xlim(0., 1.)
        axs1.set_xlabel("opinion value", fontsize=8)
        axs1.set_ylabel("number of users", fontsize=8)

        #plot media opinion evolution ..................................................................................................................................................
        axs2 = fig.add_subplot(gs[2,0])
        axs2.set_title("Opinion evolution", fontsize=11, horizontalalignment='left', x=0)
        axs2.set_ylabel("iteration step", fontsize=8)
        axs2.set_ylim(num_steps/write_every, 0)
        axs2.xaxis.grid(lw=0.5)
        axs2.set_xlim(0., 1.)
        axs2.plot(data_m[:,:], time, lw=0.01, alpha=1, color=colors['media opinion'])

        #plot final opinion state .............................................................................................
        axs3 = fig.add_subplot(gs[3,0])
        axs3.set_title("Final opinion distribution", fontsize=11, horizontalalignment='left', x=0)
        axs3.scatter(data_m[-1,:], user_numbers[-1, :], s=50, color=colors['media opinion'])
        axs3.xaxis.grid(lw=0.5)
        axs3.yaxis.grid(lw=0.5)
        axs3.set_xlabel("opinion value", fontsize=8)
        axs3.set_ylabel("number of users", fontsize=8)
        axs3.set_xlim(0., 1.)

        plt.savefig(out_path)
