import os
import logging
import warnings
from typing import Union, Dict, Callable

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

from utopya import DataManager, UniverseGroup
from utopya.plotting import UniversePlotCreator, PlotHelper, is_plot_func

# Get a logger
log = logging.getLogger(__name__)

# Increase log threshold for animation plotting

@is_plot_func(creator_type=UniversePlotCreator,
              supports_animation=False, helper_defaults=dict(
                set_title=dict(title="Opinion clusters"),
                set_labels=dict(x="Opinion", y="Step")
              ))
def densities(dm: DataManager, *,
               uni: UniverseGroup,
               hlpr: PlotHelper,
               val_range: tuple=(0., 1.),
               to_plot: str,
               num_bins: int=100):


    data        = uni['data/OpDyn/nw_users/'+to_plot]
    life_cycle  = int(uni['cfg']['OpDyn']['life_cycle'])
    time_steps  = data['time'].size
    limit   = data['time'].data[-1]
    if uni['cfg']['OpDyn']['user_ageing']=='on':
        limit   = int(limit/life_cycle)
    start, stop = val_range if val_range else (0., 1.)
    bins = num_bins
    
    data_to_plot = np.zeros((time_steps, bins))
    for row in range(time_steps):
        counts_at_time, bin_edges = np.histogram(data[row, :],
                                                 range=(start, stop),
                                                 bins=bins)
        data_to_plot[row, :] = counts_at_time/np.max(counts_at_time)

    hlpr.ax.set_xticks([i for i in np.linspace(0, bins-1, 11)])
    hlpr.ax.set_xticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
    hlpr.ax.set_yticks([i for i in np.linspace(0, time_steps-1, limit+1)])
    hlpr.ax.set_yticklabels([i for i in np.linspace(0, limit, limit+1)])
    hlpr.ax.imshow(data_to_plot[:, :], cmap='BuGn', aspect=num_bins/time_steps)
