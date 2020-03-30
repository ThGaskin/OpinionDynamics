import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from scipy.optimize import curve_fit
from utopya import DataManager, UniverseGroup
from utopya.plotting import is_plot_func, UniversePlotCreator
from .distribution_info import get_dist_info

# -----------------------------------------------------------------------------

#@is_plot_func(creator_type=UniversePlotCreator)
def user_parameters_evolution(dm: DataManager, *, out_path: str,
    uni: UniverseGroup):
        cfg         = uni['cfg']
        life_cycle  = cfg['OpDyn']['life_cycle']
        num_media   = cfg['OpDyn']['nw_m']['num_vertices']
        num_users   = cfg['OpDyn']['nw_u']['num_vertices']
        #datasets...............................................................
        age            = uni['data/OpDyn/nw_users/age_u']
        opinion        = uni['data/OpDyn/nw_users/opinion_u']
        susceptibility = uni['data/OpDyn/nw_users/susceptibility_u']
        tolerance      = uni['data/OpDyn/nw_users/tolerance_u']

        time_steps     = opinion['time'].size
        time           = opinion['time'].data
        #figure layout..........................................................
        fig     = plt.figure(figsize=(30, 40))
        widths  = [4, 4, 4]
        heights = [1.5, 0.5, 1, 2, 2, 2, 2, 2, 0.5, 1, 2, 2, 2, 2, 2]
        gs      = fig.add_gridspec(ncols=len(widths),
                                   nrows=len(heights),
                                   width_ratios=widths,
                                   height_ratios=heights,
                                   left=0.05,right=0.95,top=0.98,bottom=0.05,
                                   hspace=0)
        #colors and fonts.......................................................
        colors    = {'opinion': 'teal',
                     'tolerance':'mediumaquamarine',
                     'age': 'tomato',
                     'susceptibility':'tomato'}

        fontsizes = {'title': 60,
                     'subtitle': 40,
                     'text':30,
                     'label title': 18}

        #plots .................................................................
        for i in range(0, 5):
            j = int(i*(time_steps-1)/4)
            #opinions
            op_axs = fig.add_subplot(gs[i+3, 0])
            op_axs.hist(opinion[j,:],
                        bins=100,
                        color=colors['opinion'],
                        orientation=u'vertical')
            #tolerance
            tol_axs = fig.add_subplot(gs[i+3, 1])
            tol_axs.hist(tolerance[j,:],
                         bins=100,
                         color=colors['tolerance'],
                         orientation=u'vertical')
            #age distribution
            age_axs = fig.add_subplot(gs[i+10, 0])
            age_axs.hist(age[j,:],
                         bins=50,
                         color=colors['age'],
                         orientation=u'vertical')
            #susceptibility
            sus_axs = fig.add_subplot(gs[i+10, 1])
            sus_axs.hist(susceptibility[j,:],
                         bins=100,
                         color=colors['susceptibility'],
                         orientation=u'vertical')
            #add timestamp
            for axis in [op_axs, tol_axs, age_axs, sus_axs]:
                axis.xaxis.grid(lw=0.5)
                if not axis == age_axs:
                    axis.set_xlim((0., 1.))
                if cfg['OpDyn']['user_ageing']:
                    timestamp = str(time[j]/life_cycle)+" yrs"
                else:
                    timestamp = "step "+str(time[j])

                axis.text(0.01, 0.95, timestamp, transform=axis.transAxes,
                          fontsize=fontsizes['label title'],
                          verticalalignment='top')

        #text boxes.............................................................
        #title box
        title = fig.add_subplot(gs[0, 0:3])
        line  = lines.Line2D([0, 1], [1, 1], lw=20, color='black', axes=title)
        title.axis('off')
        title.add_line(line)
        title.text(0, 0.8, "User parameters", fontweight='bold',
                                              fontsize=fontsizes['title'],
                                              verticalalignment='top',
                                              horizontalalignment = 'left')
        subtitle = str(time[-1])+" numerical steps"
        if cfg['OpDyn']['user_ageing']:
            subtitle+="/"+str(time[-1]/life_cycle)+" years"
        subtitle+="; "+str(num_users)+" users"
        if cfg['OpDyn']['media_status']=='on':
            subtitle+="; "+str(num_media)+" media"

        title.text(0, 0.05, subtitle, fontsize=fontsizes['subtitle'])

        #subtitle boxes
        op_title       = fig.add_subplot(gs[2, 0])
        tol_title      = fig.add_subplot(gs[2, 1])
        dist_title     = fig.add_subplot(gs[2, 2])
        age_title      = fig.add_subplot(gs[9, 0])
        sus_title      = fig.add_subplot(gs[9, 1])
        subtitle_boxes = [op_title, tol_title, age_title, sus_title, dist_title]
        subtitles      = ['opinion',
                          'tolerance',
                          'age',
                          'susceptibility',
                          'distribution parameters']
        for i in range(5):
            subtitle_boxes[i].axis('off')
            line = lines.Line2D([0, 1], [0.7, 0.7], lw=7, color='black',
                                axes=subtitle_boxes[i])
            subtitle_boxes[i].text(0, 0.3, subtitles[i],
                                   fontsize=fontsizes['subtitle'],
                                   fontweight='bold')
            subtitle_boxes[i].add_line(line)

        #susceptibility distribution function
        if (get_dist_info(cfg, 'users', 'susceptibility')=='Age-dependent'):
            susc_cfg = cfg['OpDyn']['susceptibility']['users']['custom']
            s_0 = float(susc_cfg["peak"])
            s_1 = float(susc_cfg["val_at_0"])
            s_2 = float(susc_cfg["val_at_peak"])
            c = 1./s_2
            b = (1.-c*s_1)/(c*s_1*(s_0**2))
            def sus_func(arg):
                denom=c*(1.+b*((arg-s_0)**2))
                return 1/denom
            x = np.linspace(0, 100)

        #distribution info boxes
        dist_info_1 = fig.add_subplot(gs[4, 2])
        dist_info_2 = fig.add_subplot(gs[10:12, 2])
        dist_info_1.axis('off')
        dist_info_2.axis('off')
        info_1 = "Initial opinion distribution: \n"+\
                  get_dist_info(cfg, 'users', 'opinion')+"\n \n"+\
                  "Initial tolerance distribution: \n"+\
                  get_dist_info(cfg, 'users', 'tolerance')+"\n \n"+\
                  r"radicalisation parameter $k=$"+\
                  str(cfg['OpDyn']['radicalisation_parameter'])

        info_2 = "Initial age distribution: \n"+\
                  "Uniform over interval [1, 85]"+"\n \n"+\
                  "Initial susceptibility distribution: \n"+\
                  get_dist_info(cfg, 'users', 'susceptibility')+"\n \n"+\
                  "Distribution function:"+"\n"+\
                  r"$\mu(x)=\dfrac{1}{c(1+b(x-a)^2)}$"+"\n \n"+\
                  r"with $a=$"+str(s_0)+r", $b\sim$"+str(np.around(b, 5))+r", $c=$"+str(c)

        dist_info_1.text(0., 0.92, info_1, fontsize=fontsizes['text'])
        dist_info_2.text(0., 0.15, info_2, fontsize=fontsizes['text'])

        susc_plot=fig.add_subplot(gs[12:14, 2])
        susc_plot.plot(x, sus_func(x), color=colors['susceptibility'])
        susc_plot.grid(lw=1, color='slategray')

        plt.savefig(out_path)
