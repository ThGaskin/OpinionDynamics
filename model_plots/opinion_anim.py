"""This module provides plotting functions to visualize agent based models"""

import os
import logging
import warnings
from typing import Union, Dict, Callable

import numpy as np

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
                set_title=dict(title="Opinion distribution"), set_labels=dict(x="Values", y="Counts")
              ))
def opinion_animation(dm: DataManager, *,
               uni: UniverseGroup,
               hlpr: PlotHelper,
               time_idx: int,
               to_plot: str,
               num_bins: int=100,
               val_range: tuple=(0., 1.)):

    opinions    = uni['data/OpDyn/nw_users/'+to_plot]
    life_cycle  = int(uni['cfg']['OpDyn']['life_cycle'])
    time        = opinions['time'].data
    time_steps  = time.size
    if not time_idx:
        time_idx = int(next(iter([i for i in range(time_steps)])))
    if time_idx == -1:
        opinions_at_time = opinions[-1, :]

    else:
        opinions_at_time = opinions[time_idx, :]

    bins=num_bins if num_bins else 100
    start, stop = val_range if val_range else (0., 1.)

    # Calculate the opinion histogram at time_idx
    counts_at_time, bin_edges = np.histogram(opinions_at_time,
                                             range=(start, stop),
                                             bins=bins)

    # Calculate bin positions, i.e. midpoint of bin edges
    bin_pos = bin_edges[:-1] + (np.diff(bin_edges) / 2.)
    bar = hlpr.ax.bar(bin_pos, counts_at_time, width=np.diff(bin_edges),
                      color='dodgerblue')
    time_text = hlpr.ax.text(0.02, 0.95, '', transform=hlpr.ax.transAxes)
    def update_data(stepsize: int=1):
        """Updates the data of the imshow objects"""
        log.info("Plotting animation with %d frames ...",
                    opinions.shape[0] // stepsize)

        next_frame_idx = 0

        if time_steps < stepsize:
            warnings.warn("Stepsize is greater than number of steps. "
                          "Continue by plotting fist and last frame.")

        for i in range(time_steps):
            time_idx = i
            if time_idx < next_frame_idx and time_idx < time_steps:
                continue

            log.debug("Plotting frame for time index %d ...", time_idx)

            # Get the opinions
            opinion = opinions[time_idx, :]

            counts_step, _ = np.histogram(opinion, bin_edges)

            for idx, rect in enumerate(bar):
                rect.set_height(counts_step[idx])

            if(uni['cfg']['OpDyn']['user_ageing'] == 'on'):
                time_text.set_text('%.0f years' % int(time[time_idx]/life_cycle))
            else:
                time_text.set_text('step %.0f' % time[time_idx])

            hlpr.ax.relim()
            hlpr.ax.autoscale_view(scalex=False)

            next_frame_idx = time_idx + stepsize

            yield

    # Register this update method with the helper, which takes care of the rest
    hlpr.register_animation_update(update_data)
