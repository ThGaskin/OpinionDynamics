import os
import logging
import warnings
from typing import Union, Dict, Callable

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation
from matplotlib.colors import ListedColormap

from utopya import DataManager, UniverseGroup
from utopya.plotting import UniversePlotCreator, PlotHelper, is_plot_func

# Get a logger
log = logging.getLogger(__name__)

# Increase log threshold for animation plotting
logging.getLogger('matplotlib.animation').setLevel(logging.WARNING)

@is_plot_func(creator_type=UniversePlotCreator,
              supports_animation=True, helper_defaults=dict(
                set_title=dict(title="Opinion distribution by age group"),
                set_labels=dict(x="Opinion", y="Counts")
              ))
def age_groups(dm: DataManager, *,
               uni: UniverseGroup,
               hlpr: PlotHelper,
               time_idx: int,
               to_plot: str,
               num_bins: int=50,
               val_range: tuple=(0., 1.),
               ages: list):

    data        = uni['data/OpDyn/nw_users/'+to_plot]
    user_ages   = uni['data/OpDyn/nw_users/age_u']
    life_cycle  = int(uni['cfg']['OpDyn']['life_cycle'])
    time_steps  = data.coords['time'].size
    time        = data.coords['time'].data
    labels      = [(str(ages[i])+"–"+str(ages[i+1])) for i in range(len(ages)-1)]
    if (ages[-1] > 120):
        labels[-1]=str(ages[-2])+"+"
  
    if not time_idx:
        time_idx = int(next(iter(np.linspace(0, time_steps-1, 1))))

    bins = num_bins if num_bins else 50
    start, stop = val_range if val_range else (0., 1.)

    # group data by age bins
    data_by_age = np.zeros((time_steps, bins, len(ages)-1))
    for t in range (time_steps):
        data_bins=pd.cut(np.asarray(data[t, :]), bins, labels=False, include_lowest=True)
        age_bins=pd.cut(np.asarray(user_ages[t, :]), ages, labels=False, include_lowest=True)
        for i in range(len(data_bins)):
            data_by_age[t, int(data_bins[i]), int(age_bins[i])]+=1

    X=[pd.DataFrame(data_by_age[i, :, :], \
                    index=[np.around(np.linspace(start, stop, bins), 3)], \
                    columns=labels) \
                    for i in range(time_steps)]
    log.info("Created %d dataframes for plotting", len(X))

    #colors to use
    colors=['gold',
            'orange',
            'darkorange',
            'orangered',
            'firebrick',
            'darkred',
            'indigo',
            'navy',
            'royalblue',
            'cornflowerblue',
            'slategray',
            'peru',
            'saddlebrown',
            'black' ]

    if time_idx == -1:
        im = X[-1].plot.bar(stacked=True, ax=hlpr.ax, color=colors[0:-1:np.maximum(1, int(len(colors)/len(ages)))], legend=False)

    else:
        im = X[time_idx].plot.bar(stacked=True, ax=hlpr.ax, color=colors, legend=False)

    def update_data(stepsize: int=1):
        """Updates the data of the imshow objects"""
        log.info("Plotting animation with %d frames ...",
                    data.shape[0] // stepsize)

        next_frame_idx = 0

        if time_steps < stepsize:
            warnings.warn("Stepsize is greater than number of steps. "
                          "Continue by plotting fist and last frame.")
            stepsize=time_steps-1

        for time_idx in range(time_steps):
            if time_idx < next_frame_idx and time_idx < time_steps:
                continue
            hlpr.ax.clear()
            im = X[time_idx].plot.bar(stacked=True, ax=hlpr.ax, legend=False, color=colors[0:-1:np.maximum(1, int(len(colors)/len(ages)))], rot=0)
            time_text = hlpr.ax.text(0.02, 0.95, '', transform=hlpr.ax.transAxes)
            timestamp = int(time[time_idx]/life_cycle)
            time_text.set_text('%.0f years' % timestamp)
            hlpr.ax.legend(loc='upper right')
            hlpr.ax.set_xticks([i for i in np.linspace(0, bins-1, 11)])
            hlpr.ax.set_xticklabels([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
            hlpr.ax.set_title(hlpr.axis_cfg['set_title']['title'])
            hlpr.ax.set_xlabel(hlpr.axis_cfg['set_labels']['x'])
            hlpr.ax.set_ylabel(hlpr.axis_cfg['set_labels']['y'])
            next_frame_idx = time_idx + stepsize
            yield


    # Register this update method with the helper, which takes care of the rest
    hlpr.register_animation_update(update_data)
