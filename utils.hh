#ifndef UTOPIA_MODELS_OPDYN_UTILS
#define UTOPIA_MODELS_OPDYN_UTILS

#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <spdlog/spdlog.h>


namespace Utopia::Models::OpDyn::utils{

// RANDOM DISTRIBUTION UTILITY FUNCTIONS .......................................
// generate random int in range [a, b-1]
template<typename RNGType>
int get_rand_int(int a, int b, RNGType& rng) {

    std::uniform_int_distribution<int> distribution(a,b-1);
    return (int)distribution(rng);
}

// generate random double in range [a, b]
template<typename RNGType>
double get_rand_double(double a, double b, RNGType& rng) {

    std::uniform_real_distribution<double> distribution(a,b);
    return (double)distribution(rng);
}

//generate Gauss distributed random number
template<typename RNGType>
double get_Gaussian_double(double mu, double sigma, RNGType& rng) {

    std::normal_distribution<double> distribution(mu, sigma);
    return (double)distribution(rng);
}

// get random neighbour of vertex v (for directed and undirected graphs)
template<typename NWType, typename VertexDescType, typename RNGType>
auto get_rand_nb(NWType& nw, VertexDescType& v, RNGType& rng) {

    int nb_shift = get_rand_int<RNGType>(0, out_degree(v, nw), rng);
    auto nb = adjacent_vertices(v, nw).first;
    std::advance(nb, nb_shift);
    return *nb;
}
template<typename Config>
double susceptibility(const Config cfg, const int age) {
    /*! This function returns the susceptibility of a user at a given age. */
    const double s_0 = get_as<double>("peak", cfg);
    const double s_1 = get_as<double>("val_at_0", cfg);
    const double s_2 = get_as<double>("val_at_peak", cfg);
    const double c = 1./s_2;
    const double b = (1.-c*s_1)/(c*s_1*pow(s_0, 2));
    const double denom = c*(1.+b*pow((age-s_0),2));
    return 1./denom;
}

// Helper functions ............................................................
// Initialise constant distributions
template<typename ValType, typename RNGType>
ValType set_init_uniform( std::pair<ValType, ValType> interval,
                          RNGType& rng) {

    const auto log = spdlog::get("root.OpDyn");

    if (interval.first == interval.second) {
        return interval.first;
    }
    else if (interval.second > interval.first) {
        return get_rand_double(interval.first,
          interval.second, rng);
    }
    else {
        log->error("upper limit has to be higher than the lower");
        return 0;
    }
}


template<typename RNGType>
double set_init_Gauss( std::pair<double, double> distribution_info,
                       RNGType& rng) {
    /*! This function initialises a parameter with a normally distributed value in the interval [0, 1]*/
    double parameter = get_Gaussian_double(distribution_info.first,
                                    distribution_info.second, rng);
    while (parameter<=0 || parameter>=1) {
      parameter = get_Gaussian_double(distribution_info.first,
        distribution_info.second, rng);}

    return parameter;
}

// initialize properties
template <typename RNGType, typename Config>
double initialize ( const Config& cfg,
                    RNGType& rng) {

    const auto log = spdlog::get("root.OpDyn");

    std::string distribution_type =
              get_as<std::string>("distribution_type", cfg);
    if (distribution_type == "constant") {
        return get_as<double>("const_val", cfg);
    }

    if (distribution_type == "uniform") {
        std::pair<double, double> interval =
                            get_as<std::pair<double, double>>("uniform_int", cfg);
        return set_init_uniform(interval, rng);
    }

    if (distribution_type == "gaussian") {
        std::pair<double, double> distribution_info =
                           std::make_pair(get_as<double>("mean", cfg),
                                          get_as<double>("stddev", cfg));
        return set_init_Gauss(distribution_info, rng);
    }
    else {
        log->error("Invalid distribution type");
        return 0;
    }
}

/* function overloading, so the age parameter doesn't always need to be given as
an arg */
template <typename RNGType, typename Config>
double initialize (const int age,
                   const Config& cfg,
                   RNGType& rng) {

    std::string distribution_type =
         get_as<std::string>("distribution_type", cfg);

    if (distribution_type == "age-dependent") {
          const double s_0 = get_as<double>("peak", cfg["custom"]);
          const double s_1 = get_as<double>("val_at_0", cfg["custom"]);
          const double s_2 = get_as<double>("val_at_peak", cfg["custom"]);
          if (s_0<0){
            throw std::invalid_argument("Invalid value for 'peak':"
                                        "age value must be greater than 0!");
          }
          if (s_1<0 or s_1>1){
            throw std::invalid_argument("Invalid suscepbility value val_at_0:"
                                        "susceptibility must be in [0, 1]!");
          }
          if (s_2<0 or s_2>1){
            throw std::invalid_argument("Invalid suscepbility value val_at_peak:"
                                        "susceptibility must be in [0, 1]!");
          }

           return susceptibility(cfg["custom"], age);
      }

    else { return initialize(cfg, rng);}

}

} // namespace

#endif // UTOPIA_MODELS_OPDYN_UTILS
