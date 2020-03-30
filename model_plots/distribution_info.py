import numpy as np

#quickly get plot titles
def get_dist_info(config, network, parameter):
    cfg = config['OpDyn'][parameter][network]
    distribution_type = str(cfg['distribution_type'])
    if distribution_type == "constant":
        distribution_info = "Constant at " + str(cfg["const_val"])
    if distribution_type == "uniform":
        distribution_info = "Uniform over interval " + str(cfg["uniform_int"])
    if distribution_type == "gaussian":
        distribution_info = r"Gaussian ($\mu=$" + str(cfg["mean"])+r", $\sigma=$"+str(cfg["stddev"])+")"
    if distribution_type == "age-dependent":
        distribution_info = "Age-dependent"
    return distribution_info
