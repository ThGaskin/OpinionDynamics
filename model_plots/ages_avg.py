import os
import logging
import warnings
from typing import Union, Dict, Callable

import numpy as np
import pandas as pd

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
                set_title=dict(title="title"),
                set_labels=dict(x="Values", y="Counts")
              ))
def ages_avg(dm: DataManager, *,
               uni: UniverseGroup,
               hlpr: PlotHelper,
               to_plot: str,
               val_range: list=None,
               ages: list):

    data        = np.asarray(uni['data/OpDyn/nw_users/'+to_plot])
    user_ages   = uni['data/OpDyn/nw_users/age_u']
    life_cycle  = int(uni['cfg']['OpDyn']['life_cycle'])
    time_steps  = np.asarray(user_ages['time'].size)
    time        = np.asarray(user_ages['time'].data)
    user_ages   = np.asarray(user_ages)
    if uni['cfg']['OpDyn']['user_ageing']=='on':
        time=np.divide(time, life_cycle)

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

    start, stop = val_range if val_range else (0., 1.)

    hlpr.ax.set_xlim(start, stop)
    hlpr.ax.set_ylim(time[-1], 0)

    to_plot = np.zeros((time_steps, len(ages)-1))
    for t in range (0, time_steps):
        age_bins=pd.cut(user_ages[t, :], ages, labels=False, include_lowest=True)

        total_users_per_age_bin = np.zeros(len(ages)-1)
        for age in age_bins:
            total_users_per_age_bin[age]+=1

        for i in range(data.shape[1]):
            to_plot[t, age_bins[user_ages[t, i]]]+=\
                   data[t, i]/total_users_per_age_bin[age_bins[user_ages[t, i]]]

    labels = [("Ages "+str(ages[i])+" - "+str(ages[i+1])) for i in range(len(ages)-1)]
    labels[-1]="Ages "+str(ages[-2])+"+"

    for i in range(len(ages)-1):
        hlpr.ax.plot(to_plot[:, i], time,
                     lw=1, label=labels[i],
                     color=colors[0:-1:np.maximum(1, int(len(colors)/len(ages)))][i])
    hlpr.ax.legend(loc=1)
