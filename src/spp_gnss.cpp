#include "spp_gnss.h"
#include "spp.h"
void spp_gnss_process(const std::vector<std::vector<ObsData>>& observations,
                      const std::vector<NavData>& navigation, const double* reference) {
    processObservations(observations,navigation,reference,"GCER");
}
