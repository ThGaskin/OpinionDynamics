#ifndef UTOPIA_MODELS_OPDYN_UPDATE
#define UTOPIA_MODELS_OPDYN_UPDATE

namespace Utopia::Models::OpDyn::update{

// UPDATE UTILITY FUNCTIONS ....................................................
//user-user opinion update
template <typename VertexDescType, typename NWType>
void opinion( VertexDescType& v,
                     VertexDescType& nb,
                     NWType& nw)
{
    if (fabs(nw[v].opinion-nw[nb].opinion) <= nw[v].tolerance) {
        nw[v].opinion += nw[v].susceptibility * (nw[nb].opinion-nw[v].opinion);
    }
}

//user-media opinion update
template <typename VertexDescType, typename NWType_1, typename NWType_2>
void opinion( VertexDescType& v,
                     VertexDescType& nb,
                     NWType_1& nw_1,
                     NWType_2& nw_2)
{
    if (fabs(nw_1[v].opinion-nw_2[nb].opinion) <= nw_1[v].tolerance) {
        nw_1[v].opinion += nw_1[v].susceptibility * nw_2[nb].persuasiveness
                            * (nw_2[nb].opinion-nw_1[v].opinion);
    }
}

template <typename VertexDescType, typename NWType>
void tolerance(VertexDescType& v,
                      NWType& nw,
                      const double previous_opinion,
                      const double radicalisation_parameter)
{
    const double x = fabs(previous_opinion-0.5);
    const double y = fabs(nw[v].opinion-0.5);
    const double r = 1.+radicalisation_parameter*(pow(y, 2)-pow(x, 2));
    nw[v].tolerance=pow(nw[v].tolerance, r);
}

} // namespace

#endif // UTOPIA_MODELS_OPDYN_UPDATE
