## Fundamentals

Opinion dynamics aim to simulate and analyse the spread of opinions on adaptive networks. Each node in the network represents a user, and each user is connected through edges to neighbours with whom the user can interact. This specific version of an opinion dynamics model is based on Deffuant and Weisbruch's [Deffuant 2000] selective exposure model: users do not interact with the entirety of the opinion space, but only with users whose opinions fall within the user's **tolerance range**.

In this model, we select a single random user in each time step, and perform an opinion update, that is, we let the user interact with one of her neighbours. The neighbour is chosen according to the weight distribution of each edge, which are normalised to represent an interaction probability. An edge weight of 0.5 thus represents a 50% chance of that neighbour being selected for interaction in a time step. The user network is a **directed network**, meaning in- and out-edges have different weights.

In addition to the user network, we have implemented a **media network** (see eg. [Quattrociocchi 2014]) to simulate the interaction of people with influencers (politicians, media companies, celebrities, etc.) -  users or entities, whose sole purpose is to attract as many followers as possible, and who will adjust their opinion to maximise their user count.

## Implementation

### Users
Each user u<sub>i</sub> is given a set of five parameters: **opinion &sigma;<sub>i</sub>**, **tolerance &epsilon;<sub>i</sub>**, **susceptibility &mu;<sub>i</sub>**, **age a**, and (if the media network is turned on) the **used medium**.

### User opinion
The opinion lives in the **opinion space**, which in this model is the unit interval [0, 1] as a subset of R. The opinion values are thus real numbers (doubles). The initial opinion distribution can be intialised as a constant, uniformly distributed, Gauss-distributed, or age-distributed (cf. the <code>distribution_type</code> key below). They are updated according to the selective exposure law

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_i(t&plus;1)=\sigma_i(t)&plus;\mu_i&space;*&space;(\sigma_j(t)-\sigma_i(t))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_i(t&plus;1)=\sigma_i(t)&plus;\mu_i&space;*&space;(\sigma_j(t)-\sigma_i(t))" title="\sigma_i(t+1)=\sigma_i(t)+\mu_i * (\sigma_j(t)-\sigma_i(t))" /></a>     &nbsp;&nbsp;&nbsp; (1)

if

<a href="https://www.codecogs.com/eqnedit.php?latex=|\sigma_j(t)-\sigma_i(t)|&space;\leq&space;\epsilon_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?|\sigma_j(t)-\sigma_i(t)|&space;\leq&space;\epsilon_i" title="|\sigma_j(t)-\sigma_i(t)| \leq \epsilon_i" /></a>.

### User tolerance
The user tolerance is a double value in [0, 1]. It decreases when users radicalise (ie. move away from the "moderate opinion" 0.5) and increases when they deradicalise. The radicalisation law is

<a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon_i&space;(t&plus;1)&space;=&space;\epsilon_i&space;(t)^{f(\sigma_i&space;(t&plus;1),&space;\&space;\sigma_i(t))}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon_i&space;(t&plus;1)&space;=&space;\epsilon_i&space;(t)^{f(\sigma_i&space;(t&plus;1),&space;\&space;\sigma_i(t))}" title="\epsilon_i (t+1) = \epsilon_i (t)^{f(\sigma_i (t+1), \ \sigma_i(t))}" /></a>  &nbsp;&nbsp;&nbsp; (2),


<a href="https://www.codecogs.com/eqnedit.php?latex=f(x,&space;y)=1-k((x-0.5)^2-(y-0.5)^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f(x,&space;y)=1-k((x-0.5)^2-(y-0.5)^2)" title="f(x, y)=1-k((x-0.5)^2-(y-0.5)^2)" /></a>  &nbsp;&nbsp;&nbsp; (3).

The **radicalisation parameter** k determines the strength of the effect, and must be a value in [0, 4]. Setting k=0 turns off the tolerance change effect entirely. Note that this function skews the tolerance distribution towards low values - this is an implementation error that is yet to be corrected.

### User susceptibility
A user's susceptibility parameter can be set as a function of her age a(t):

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu(x)&space;=&space;\dfrac{1}{c(1&plus;b(x-a)^2)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu(x)&space;=&space;\dfrac{1}{c(1&plus;b(x-a)^2)}" title="\mu(x) = \dfrac{1}{c(1+b(x-a)^2)}" /></a> &nbsp;&nbsp;&nbsp; (4),

where a, b, and c are parameters that specify the shape of the curve. In the configuration file, these parameters are transformed into more user-friendly variables (age of peak susceptibility, susceptibility at peak, susceptibility at age 0). This gives users complete control over the susceptibility distribution function.

### User ageing
A user's age is initialised at random. If ageing is turned on, every year (the length of year having been set in the config-file) a fraction of the senior members of the population are reinitialised as young users of age 1. They are rewired to a parent node, as well as to several peers (peer, children, and senior ages can be set in the config). The child-parent edge is given weight 0.5, the out edges to the peers are given weight 1/(out degree -1). The child's tolerance is set to the parent's tolerance, and the child's opinion value is set a mixture of the parent opinion and the average opinion of the child's peers.

### Weights
After each opinion revision, the user can cut links to a neighbour whose opinion is furthest from her own. A new neighbour is chosen at random, and the weights set to 1/(number of users added). Additionally, weights are adjusted in each time step to reflect two things: distance in opinion space, and the age difference. Weights are adjusted according to

<a href="https://www.codecogs.com/eqnedit.php?latex=w_{ij}(t&plus;1)=w_{ij}(t)&space;\cdot&space;\left(1-w&space;*&space;|\sigma_j&space;-&space;\sigma_i|&space;&plus;&space;e^{\log(0.5)&space;\cdot&space;|\Delta&space;a|/(2a)}\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?w_{ij}(t&plus;1)=w_{ij}(t)&space;\cdot&space;\left(1-w&space;*&space;|\sigma_j&space;-&space;\sigma_i|&space;&plus;&space;e^{\log(0.5)&space;\cdot&space;|\Delta&space;a|/(2a)}\right&space;)" title="w_{ij}(t+1)=w_{ij}(t) \cdot \left(1-w * |\sigma_j - \sigma_i| + e^{\log(0.5) \cdot |\Delta a|/(2a)}\right )" /></a>  &nbsp;&nbsp;&nbsp; (5),

where the **weighting parameter** determines the strength of the selective exposure and &Delta;a is the age difference between the two users in question. Users thus prefer neighbours whose age and opinion is close to their own.

### Media network
The media network is an **undirected network** between competing media nodes. Each medium also holds an opinion, has a tolerance radius, is susceptible to opinion change, and knows how many followers (voters, readers, ...) it currently has. Additionally, a medium has a **persuasiveness &rho;**, which describes how convincing a medium is to its user base, and an advertisement budget **ads**, which describes how much money a medium can spend on outreach to attract new followers. The persuasiveness determines the strength of the media-user network coupling. A user interacts with a medium via the interaction law

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_i(t&plus;1)&space;=&space;\sigma_i(t)&plus;\mu_i&space;\cdot&space;\rho_j&space;\cdot&space;(\sigma_j(t)-\sigma_i(t))" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_i(t&plus;1)&space;=&space;\sigma_i(t)&plus;\mu_i&space;\cdot&space;\rho_j&space;\cdot&space;(\sigma_j(t)-\sigma_i(t))" title="\sigma_i(t+1) = \sigma_i(t)+\mu_i \cdot \rho_j \cdot (\sigma_j(t)-\sigma_i(t))" /></a>  &nbsp;&nbsp;&nbsp; (6),

where the index j represents the medium properties. In each media update, a single medium compares its own user count to that of its neighbours, and shifts its own stance to that of the most popular competitor, but only if the opinion distance is smaller than the tolerance, and larger than third of the medium's tolerance. If the opinion distance is less than a third of the tolerance, the medium moves away from the competitor. The size of a medium's user base directly leads to ad revenue. The higher a medium's advertisement budget, the more users are exposed to that medium's opinion. In each time step, a random user is assigned a particular medium based on that medium's ad value, and performs an opinion update if the medium's opinion is within her tolerance range.
## How to run the model

1. <code>model</code>: Choose the network type. For the directed user network <code>nw_u</code>, available models are <code>ErdosRenyi</code> (random), <code>BollobasRiordan</code> (scale-free directed), or <code>WattsStrogatz</code> (small-world). For the undirected media network <code>nw_m</code>, the scale-free option requires the <code>BarabasiAlbert</code> key.
2. <code>num_vertices</code>: Select number of users or media.
3. <code>mean_degree</code>: Mean number of out-edges of a single node.
4. <code>opinion</code>, <code>tolerance</code>, <code>susceptibility</code>, <code>persuasiveness</code>. This sets the initial distributions of the various vertex properties for <code>users</code> and <code>media</code>:
* <code>distribution_type</code>: Choose from <code>constant</code>, <code>uniform</code>, <code>gaussian</code>, or <code>age-dependent</code> (only for user susceptibility).
* <code>const_val</code>: If the distribution type is constant, set the constant value here.
* <code>uniform_int</code>: If the distribution type is uniform, specify the interval over which the parameter is to be set to a uniform value.
* <code>mean</code>/<code>stddev</code>: If the distribution type is Gaussian, specify the mean and standard deviation.

For user susceptibility, a custom susceptibility function can be defined via the <code>custom</code> key:

* <code>peak</code>: Specify the age of peak susceptibility (must be greater or equal to 0!).
* <code>val_at_peak</code>: Specify the susceptibility at the peak (must in [0, 1]).
* <code>val_at_0</code>: Specify the susceptibility at age 0 (must be in [0, 1]).
5. <code>radicalisation_parameter</code>: Determines the strength of the radicalisation effect (see equation 3). Set it to 0 to turn off all radicalisation effects. Cannot be larger than 4, or else the tolerance can become greater than 1.
6. <code>media_status</code>: Turn the media network either to <code>on</code> or <code>off</code>.
7. <code>life_cycle</code>: Determines how many numerical steps are equal to a single year. If the <code>user_ageing</code> key is turned off, this key has no effect.
8. <code>replacement_rate</code>: Determines how many users are reinitialised as children each year; 2-3% are realistic values.
9. <code>age_groups</code>: Set the various age groups for <code>children</code>, <code>parents</code>, and <code>seniors</code>.
10.<code>media_time_constant</code>: If the media network is turned on, it can run on a slower timescale than the user network, ie. media will update their opinions less frequently than the public. Setting this key to a value of 2, for instance, lets the media network update itself every second time step. Must be an integer larger than 0.
11. <code>init_ads</code>: The initial ad value interval.
12. <code>weighting</code>: The weighting parameter from equation (5).
10. <code>rewiring</code>: The probability that a user will rewire ties to neighbours furthest away in opinion space.

### Output
The model outputs several user data plots and one media data plot (if the media network is turned on):

#### 1. User parameters plot:

![user_param_plot](https://ts-gitlab.iup.uni-heidelberg.de/uploads/-/system/personal_snippet/27/4319df8fd03a54560942f159e89672c7/user_parameters.png)

This plot shows user opinion, tolerance, age group distribution, and susceptibility for five timeframes. The column on the right displays the model configuration values, as well as the susceptibility distribution function (if <code>distribution_type = age-dependent</code> is selected.)

#### 2. Animated user opinion evolution:

![op_plot](https://ts-gitlab.iup.uni-heidelberg.de/uploads/-/system/personal_snippet/26/10e403b5e167c6da42682c7ae41d80ed/op_plot.png)

This animation shows the opinion distribution over time.

#### 3. Animated opinion distribution by age groups:

![age_groups](https://ts-gitlab.iup.uni-heidelberg.de/uploads/-/system/personal_snippet/34/3ae51f3288fb83ce1d4e212cb673bfd8/age_groups.png)

This animated stacked bar plot shows the opinion distribution for various age groups (the age groups to be plotted can be specified in the <code>plot_cfg.yml</code> file)
Finally, when the media network is being used, the model outputs a third plot:

##### 4. Average opinion of age groups:

![age_grp_avg](https://ts-gitlab.iup.uni-heidelberg.de/uploads/-/system/personal_snippet/35/089319002bd34bde74d520fa39976234/age_grp_avgs.png)

For each age group (can be specified in the plot configuration) the average opinion is plotted over time.

##### 5. Opinion clusters
![op_cluster](https://ts-gitlab.iup.uni-heidelberg.de/uploads/-/system/personal_snippet/36/7c7902d705c613a4bd86aeea77b6bc44/op_clust.png)

This plot outputs the opinion clusters over time: for each time step, the graph plots the opinion values divided by the maximum value in that time step.

##### 6. Media plot

![media](https://ts-gitlab.iup.uni-heidelberg.de/uploads/-/system/personal_snippet/28/0f8a3c3414f4d6a0cb91632f4be5671a/media.png)

The scatterplots show the initial media distribution (user numbers over opinion value), the middle plot displays the media evolution.

## Future projects
- Update users individually and throughout the year, rather than in batches
- Model user properties (especially age-dependent ones) more closely on real sociological/psychological empirical findings
- Generalise the concept of "radicalisation" to refer to average user behaviour rather than absolute position on opinion spectrum
- Long-term plan: opinion dynamics on different opinion space topologies!

## References
- Deffuant G. et al: _Mixing beliefs among interacting agents._ Adv Complex Syst. (2000) 3:87-98. DOI: 10.1142/S0219525900000078
- Quattrociocchi W. et al: _Opinion dynamics on interacting networks: media competition and social influence._ Sci. Rep. (2014) 4, 4938. DOI: 10.1038/srep04938
