#ifndef MYGNSS10_SPP_H
#define MYGNSS10_SPP_H
#include "rinex.h"
#include <map>
#include <string>
#include <utility>

struct SPPResult {
    GPSTime time{};
    double pos[3]{}, vel[3]{}; // ECEF m and m/s
    std::map<char,double> clockBias, clockDrift; // s and s/s, per system
    int nSat=0, nVel=0;
    int nRejected=0; // pseudorange outliers removed before the final solution
    double pdop=0;
    bool positionValid=false, velocityValid=false;
    std::string status="no_observations";
};

void ecef2pos(const double* r, double* pos);
void ecefVectorToEnu(const double* v, const double* reference, double* enu);
bool observationFrequencies(const ObsData& obs, const NavData& nav, double& f1, double& f2);

// The solver owns its navigation index; no dangling pointers to the caller's vector.
class SPPSolver {
public:
    explicit SPPSolver(const std::vector<NavData>& navigation);
    SPPResult solve(const std::vector<ObsData>& epoch, const double* initialPosition,
                    const std::string& systems="GCER") const;
    // For diagnostics: same signal/health/age selection as the positioning solver.
    // Returned pointer remains valid for this solver's lifetime.
    const NavData* selectEphemeris(const ObsData& observation) const;
private:
    std::map<std::pair<char,int>,std::vector<NavData>> navigation_;
};
void spp_process(const std::vector<std::vector<ObsData>>& observations,
                 const std::vector<NavData>& navigation, const double* reference);
int processObservations(const std::vector<std::vector<ObsData>>& observations,
                         const std::vector<NavData>& navigation, const double* reference,
                         const std::string& systems, bool coldStart=false);
#endif
