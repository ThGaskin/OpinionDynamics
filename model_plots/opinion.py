import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from utopya import DataManager, UniverseGroup
from utopya.plotting import is_plot_func, UniversePlotCreator
from .distribution_info import get_dist_info

# -----------------------------------------------------------------------------

#@is_plot_func(creator_type=UniversePlotCreator)
def user_opinion_evolution(dm: DataManager, *, out_path: str,
            uni: UniverseGroup):
    """
    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (UniverseGroup): The selected universe data
    """
    #get data
    uni_cfg         = uni['cfg']
    num_steps       = uni_cfg['num_steps']
    life_cycle      = uni_cfg['OpDyn']['life_cycle']
    write_every     = uni_cfg['write_every']
    num_users       = uni_cfg['OpDyn']['nw_u']['num_vertices']
    avg_user_degree = uni_cfg['OpDyn']['nw_u']['mean_degree']
    num_media       = uni_cfg['OpDyn']['nw_m']['num_vertices']
    data_u          = uni['data/OpDyn/nw_users/opinion_u']
    age_u           = uni['data/OpDyn/nw_users/age_u']
    time            = opinion['time'].data

    #set up format
    fig = plt.figure(figsize=(40,20), constrained_layout=True)
    widths = [6, 3]
    heights = [1, 1, 1, 1, 1, 1]
    gs = fig.add_gridspec(ncols=len(widths), nrows=len(heights), width_ratios=widths, height_ratios=heights)

    #colors and fonts
    colors = {'opinion evolution':'steelblue',
              'time plot':'dodgerblue'}

    fonts = {'title': 60,
             'subtitle': 30,
             'text':25,
             'plot title': 20,
             'label title': 18}

    alpha = np.exp(-num_users/1000.)

    #set binning and range
    num_bins = 100
    space=(0., 1.)

    # calculate histogram arrays and age distributions
    def get_hist_data(timestep):
        counts, bin_edges = np.histogram(data_u[timestep,:], num_bins, space)
        average_user_age = np.zeros(num_bins)
        x_axis = np.zeros(num_bins)
        bins = np.linspace(space[0], space[1], num_bins+1)
        indices = np.digitize(data_u[timestep,:], bins)
        for i in range(len(data_u[timestep,:])):
            average_user_age[indices[i]-1]+=age_u[timestep, i]
        for i in range(len(average_user_age)):
            if counts[i]!=0:
                average_user_age[i]/=counts[i]
            x_axis[i]=0.5*(bin_edges[i+1]+bin_edges[i])
        data=np.vstack((x_axis, counts, average_user_age))
        return data

    #assign each age bin a color if user ageing is turned on
    def get_age_color(timestep):
        if uni_cfg['OpDyn']['user_ageing']:
            res = []
            age_array = get_hist_data(timestep)[2, :]
            normalising_factor = np.max(get_hist_data(timestep)[2, :])
            for i in range(len(age_array)):
                col = age_array[i]/normalising_factor
                res.append((col, 0.2, 1.-col, 1.))
            return res
        else:
            return colors['time plot']

    #axis 0 is just for text
    axs0 = fig.add_subplot(gs[0, 0])
    axs0.axis('off')
    axs0.text(0, 1.0, "User opinion evolution", fontweight='bold', fontsize=fonts['title'], verticalalignment='top', horizontalalignment = 'left')
    if uni_cfg['OpDyn']['user_ageing']:
        axs0.text(0, 0.5, str(num_steps)+" numerical steps, corresponding to "+str(num_steps/life_cycle)+" years", fontweight='bold', fontsize=fonts['text'])
    else:
        axs0.text(0, 0.5, "Numerical steps: "+str(num_steps), fontweight='bold', fontsize=fonts['text'])
    axs0.text(0, 0.35, "Number of users: "+str(num_users), fontweight='bold', fontsize=fonts['text'])
    if uni_cfg['OpDyn']['media_status']:
        axs0.text(0., 0.2, "Number of media: "+str(num_media), fontweight='bold', fontsize=fonts['text'])
    axs0.text(0, 0.0, "Initial distribution: "+get_dist_info(uni_cfg, 'users', 'opinion'), fontsize=fonts['text'])

    #plot opinion evolution .............................................................................................
    axs1 = fig.add_subplot(gs[1:6,0])
    axs1.set_title("opinion evolution", fontsize=fonts['text'], horizontalalignment='left', x=0)
    axs1.set_xlabel("opinion value", fontsize=fonts['label title'])
    axs1.set_ylabel("iteration step", fontsize=fonts['label title'])
    axs1.set_xlim(0., 1.)
    axs1.set_ylim(num_steps/write_every, 0)
    axs1.xaxis.grid(lw=0.5)
    axs1.plot(data_u[:,:], time, lw=0.01, alpha=alpha, color=colors['opinion evolution'])

    #plot initial opinion state  .................................................................................................................................................
    axs2 = fig.add_subplot(gs[0, 1])
    axs2.xaxis.grid(lw=0.5)
    axs2.set_xlim(space)
    axs2.set_ylabel("user group size", fontsize=fonts['label title'])
    axs2.bar(get_hist_data(0)[0, :], get_hist_data(0)[1, :], (space[1]-space[0])/num_bins, color=(get_age_color(0))[:])

    #plot snapshots .............................................................................................................................................
    for i in range(1, 5):
          j = int((i)*num_steps/write_every/5)
          axs=fig.add_subplot(gs[i, 1])
          axs.set_ylabel("user group size", fontsize=fonts['label title'])
          axs.xaxis.grid(lw=0.5)
          axs.set_xlim(space)
          axs.bar(get_hist_data(j)[0, :], get_hist_data(j)[1, :], (space[1]-space[0])/num_bins, color=(get_age_color(j)[:]))


    #plot final opinion state .............................................................................................
    axs6 = fig.add_subplot(gs[5,1])
    axs6.set_xlabel("opinion value", fontsize=fonts['label title'])
    axs6.set_ylabel("user group size", fontsize=fonts['label title'])
    axs6.xaxis.grid(lw=0.5)
    axs6.set_xlim(space)
    axs6.bar(get_hist_data(-1)[0, :], get_hist_data(-1)[1, :], (space[1]-space[0])/num_bins, color=(get_age_color(-1))[:])


    plt.savefig(out_path)
