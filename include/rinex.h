#ifndef MYGNSS10_RINEX_H
#define MYGNSS10_RINEX_H
#include <limits>
#include <string>
#include <vector>
#include "gtime.h"

// All internal epochs are GPST. toeSow retains the native broadcast time scale.
struct NavData {
    int sat = 0;
    char sys = 0;
    GPSTime toe{}, toc{};
    double toeSow = 0;
    double a0 = 0, a1 = 0, a2 = 0;
    double crs = 0, delta_n = 0, M0 = 0;
    double cuc = 0, e = 0, cus = 0, sqrtA = 0;
    double cic = 0, Omega0 = 0, cis = 0;
    double i0 = 0, crc = 0, omega = 0, OmegaDot = 0, IDot = 0;
    double tgd = 0, tgd2 = 0;
    int health = 0, dataSources = 0, freq_num = 0;
    double gloPos[3]{}, gloVel[3]{}, gloAcc[3]{}; // m, m/s, m/s^2
};
struct ObsData {
    GPSTime time{};
    int sat = 0;
    char sys = 0;
    double P1 = 0, P2 = 0; // m; zero means missing
    // Zero Doppler is valid. Missing observations are represented separately.
    double D1 = std::numeric_limits<double>::quiet_NaN();
    double D2 = std::numeric_limits<double>::quiet_NaN();
    std::string codeP1, codeP2;
};
// Invalid/truncated input throws std::runtime_error.
std::vector<NavData> readNavFile(const std::string& filename);
std::vector<std::vector<ObsData>> readObsFile(const std::string& filename, double* approxPos);
#endif
