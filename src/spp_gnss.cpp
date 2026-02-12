/*
 * src/spp_gnss.cpp - 多系统双频 IF 组合 SPP (BDS TGD 修复版)
 * 修正记录:
 * 1. [CRITICAL] 针对 BDS 增加 TGD (Timing Group Delay) 改正。
 * 原因: BDS 广播钟差基于 B3，而 IF 组合使用 B1/B3，必须扣除 TGD1 * Alpha。
 * GPS 不需要此修正(基于 L1/L2 IF)，Galileo 类似 GPS。
 * 2. 保持之前的定权(sin^2)和坐标转换修复。
 */

#include "spp_gnss.h"
#include "gtime.h"
#include "matrix.h"
#include "orbit.h"
#include "tropo.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <cstdio>

using namespace std;

// =================================================================
// [调试开关]
// =================================================================
const bool USE_GPS = 1;
const bool USE_BDS = 1;
const bool USE_GAL = 1;
const bool USE_GLO = 1;
// =================================================================

const double CLIGHT_GNSS = 299792458.0;
const double PI = 3.1415926535897932;
const double OMEGA_E = 7.2921151467E-5;
const double D2R = PI / 180.0;
const double EL_MASK = 10.0 * D2R;

// 椭球参数
const double RE_WGS84 = 6378137.0;
const double FE_WGS84 = 1.0 / 298.257223563;

// 频率
const double FREQ_GPS_L1 = 1.57542E9;
const double FREQ_GPS_L2 = 1.22760E9;
const double FREQ_BDS_B1 = 1.561098E9;
const double FREQ_BDS_B3 = 1.26852E9;
const double FREQ_GAL_E1 = 1.57542E9;
const double FREQ_GAL_E5a = 1.17645E9;
const double FREQ_GLO_G1_BASE = 1.60200E9;
const double FREQ_GLO_G1_DELTA = 0.5625E6;
const double FREQ_GLO_G2_BASE = 1.24600E9;
const double FREQ_GLO_G2_DELTA = 0.4375E6;

static bool is_sys_enabled(char sys) {
    if (sys == 'G') return USE_GPS;
    if (sys == 'C') return USE_BDS;
    if (sys == 'E') return USE_GAL;
    if (sys == 'R') return USE_GLO;
    return false;
}

static bool get_sys_freqs(char sys, int k, double& f1, double& f2) {
    switch(sys) {
        case 'G': f1 = FREQ_GPS_L1; f2 = FREQ_GPS_L2; return true;
        case 'C': f1 = FREQ_BDS_B1; f2 = FREQ_BDS_B3; return true;
        case 'E': f1 = FREQ_GAL_E1; f2 = FREQ_GAL_E5a; return true;
        case 'R':
            f1 = FREQ_GLO_G1_BASE + k * FREQ_GLO_G1_DELTA;
            f2 = FREQ_GLO_G2_BASE + k * FREQ_GLO_G2_DELTA;
            return true;
        default: return false;
    }
}

// ECEF -> BLH
static void ecef2pos(const double *r, double *pos) {
    double e2 = FE_WGS84 * (2.0 - FE_WGS84);
    double r2 = r[0]*r[0] + r[1]*r[1];
    double z = r[2];
    double zk = 0.0;
    double v = RE_WGS84;
    double sinp = 0.0;
    for (double z_old = -1e16; fabs(zk - z_old) > 1e-4;) {
        z_old = zk;
        sinp = z / sqrt(r2 + z * z);
        v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);
        zk = r[2] + v * e2 * sinp;
    }
    pos[0] = (r2 > 1e-12) ? atan(zk / sqrt(r2)) : (r[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
    pos[1] = (r2 > 1e-12) ? atan2(r[1], r[0]) : 0.0;
    pos[2] = sqrt(r2 + z * z) - v;
}

static void xyz2enu(const double *r, const double *ref, double *enu) {
    double lat = atan2(ref[2], sqrt(ref[0]*ref[0]+ref[1]*ref[1]));
    double lon = atan2(ref[1], ref[0]);
    double sinP=sin(lat), cosP=cos(lat);
    double sinL=sin(lon), cosL=cos(lon);
    double dx=r[0]-ref[0], dy=r[1]-ref[1], dz=r[2]-ref[2];
    enu[0] = -sinL*dx + cosL*dy;
    enu[1] = -sinP*cosL*dx - sinP*sinL*dy + cosP*dz;
    enu[2] =  cosP*cosL*dx + cosP*sinL*dy + sinP*dz;
}

static void vel2enu(const double *v, const double *ref, double *ve) {
    double lat = atan2(ref[2], sqrt(ref[0]*ref[0]+ref[1]*ref[1]));
    double lon = atan2(ref[1], ref[0]);
    double sinP=sin(lat), cosP=cos(lat);
    double sinL=sin(lon), cosL=cos(lon);
    ve[0] = -sinL*v[0] + cosL*v[1];
    ve[1] = -sinP*cosL*v[0] - sinP*sinL*v[1] + cosP*v[2];
    ve[2] =  cosP*cosL*v[0] + cosP*sinL*v[1] + sinP*v[2];
}

static void sat_azel(const double *r, const double *rs, double *azel) {
    double enu[3]; xyz2enu(rs, r, enu);
    double d = sqrt(enu[0]*enu[0] + enu[1]*enu[1] + enu[2]*enu[2]);
    if (d < 1e-3) { azel[0]=0; azel[1]=0; return; }
    azel[1] = asin(enu[2] / d);
    azel[0] = atan2(enu[0], enu[1]);
    if (azel[0] < 0) azel[0] += 2*PI;
}

static double calc_IF(double obs1, double obs2, double f1, double f2) {
    double f1_sq = f1 * f1;
    double f2_sq = f2 * f2;
    return (f1_sq * obs1 - f2_sq * obs2) / (f1_sq - f2_sq);
}

// --------------------------------------------------------------------------
// 1. 位置解算
// --------------------------------------------------------------------------
static int estPosition(const vector<ObsData>& obsList, const vector<NavData>& navList,
                       double* rr, map<char, double>& clk, double* pdop) {

    vector<char> syss;
    for(const auto& o:obsList) {
        if(!is_sys_enabled(o.sys)) continue;
        if(o.P1!=0 && o.P2!=0) {
            bool found=false; for(char s:syss) if(s==o.sys) found=true;
            if(!found) syss.push_back(o.sys);
        }
    }
    if(syss.empty()) return 0;

    int m = 3 + syss.size();
    for(char s:syss) if(clk.find(s)==clk.end()) clk[s]=0.0;

    int nSat = 0;
    for(int iter=0; iter<10; iter++) {
        vector<double> H, v, W;
        nSat = 0;

        bool valid_pos = (sqrt(rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]) > 1.0);

        double pos_blh[3] = {0};
        if(valid_pos) ecef2pos(rr, pos_blh);

        for(const auto& obs : obsList) {
            if(!is_sys_enabled(obs.sys)) continue;
            if(obs.P1==0 || obs.P2==0) continue;

            const NavData* nav = nullptr; double min_dt=1e9;
            for(const auto& n:navList) {
                if(n.sat==obs.sat && n.sys==obs.sys) {
                    double dt = timediff(obs.time, n.toe);
                    if(fabs(dt)<fabs(min_dt)) { min_dt=dt; nav=&n; }
                }
            }
            if(!nav || fabs(min_dt)>14400) continue;

            double f1, f2;
            if(!get_sys_freqs(obs.sys, nav->freq_num, f1, f2)) continue;

            double tau = obs.P1 / CLIGHT_GNSS;
            GPSTime t_tx = obs.time; t_tx.sec -= tau;

            double rs[3]={0}, dts=0.0;
            satPos(t_tx, *nav, rs, &dts);

            // =======================================================
            // [FIX] BDS TGD 改正 (关键)
            // GPS L1/L2 IF 不需要 (广播钟差已包含)
            // BDS B1/B3 IF 需要 (广播钟差基于 B3)
            // =======================================================
            if (obs.sys == 'C') {
                double gamma = (f1 * f1) / (f1 * f1 - f2 * f2);
                dts -= nav->tgd * gamma;
            }
            // =======================================================

            double theta = OMEGA_E * tau;
            double rx = rs[0]*cos(theta) + rs[1]*sin(theta);
            double ry = -rs[0]*sin(theta) + rs[1]*cos(theta);
            rs[0]=rx; rs[1]=ry;

            double r = sqrt(pow(rs[0]-rr[0],2)+pow(rs[1]-rr[1],2)+pow(rs[2]-rr[2],2));

            double azel[2]={0};
            if(valid_pos) {
                sat_azel(rr, rs, azel);
                if(azel[1] < EL_MASK) continue;
            } else {
                azel[1] = PI/2.0;
            }

            double trop=0, var_trop=0;
            if(iter>0 && valid_pos) {
                trop_model_prec(obs.time, pos_blh, azel, &trop, &var_trop);
            }

            double P_IF = calc_IF(obs.P1, obs.P2, f1, f2);
            double res = P_IF - (r + clk[obs.sys] - CLIGHT_GNSS*dts + trop);

            H.push_back((rr[0]-rs[0])/r);
            H.push_back((rr[1]-rs[1])/r);
            H.push_back((rr[2]-rs[2])/r);
            for(char s:syss) H.push_back(s==obs.sys?1.0:0.0);

            v.push_back(res);

            // 定权 (sin^2)
            double w = 1.0;
            if (valid_pos) {
                double sin_el = sin(azel[1]);
                w = sin_el * sin_el;
                if (obs.sys == 'C' && obs.sat <= 5) w *= 0.1; // GEO 降权
            }
            W.push_back(w);
            nSat++;
        }

        if(nSat < m) return 0;

        vector<double> Q(m*m,0), b(m,0), Q_geo(m*m,0);
        for(int i=0;i<nSat;i++) {
            for(int j=0;j<m;j++) {
                for(int k=0;k<m;k++) {
                    Q[j*m+k] += H[i*m+j]*W[i]*H[i*m+k];
                    Q_geo[j*m+k] += H[i*m+j]*H[i*m+k];
                }
                b[j] += H[i*m+j]*W[i]*v[i];
            }
        }
        if(matinv(Q.data(), m)==-1) return 0;
        vector<double> dx(m);
        matmul("NN", m, m, 1, 1.0, Q.data(), b.data(), 0.0, dx.data());

        rr[0]+=dx[0]; rr[1]+=dx[1]; rr[2]+=dx[2];
        for(int k=0;k<(int)syss.size();k++) clk[syss[k]]+=dx[3+k];

        if(matinv(Q_geo.data(), m)!=-1) {
            *pdop = sqrt(Q_geo[0] + Q_geo[m+1] + Q_geo[2*m+2]);
        }

        if(sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2])<1e-4) return nSat;
    }
    return nSat;
}

// --------------------------------------------------------------------------
// 2. 速度解算
// --------------------------------------------------------------------------
static int estVelocity(const vector<ObsData>& obsList, const vector<NavData>& navList,
                       const double* rr, double* vr, map<char, double>& clk_drift) {

    vector<char> syss;
    for(const auto& o:obsList) {
        if(!is_sys_enabled(o.sys)) continue;
        if(o.D1!=0 && o.D2!=0) {
            bool found=false; for(char s:syss) if(s==o.sys) found=true;
            if(!found) syss.push_back(o.sys);
        }
    }
    if(syss.empty()) return 0;

    int m = 3 + syss.size();
    vr[0]=0; vr[1]=0; vr[2]=0;
    for(char s:syss) if(clk_drift.find(s)==clk_drift.end()) clk_drift[s]=0.0;

    int nSat = 0;
    for(int iter=0; iter<10; iter++) {
        vector<double> H, v, W;
        nSat = 0;

        for(const auto& obs : obsList) {
            if(!is_sys_enabled(obs.sys)) continue;
            if(obs.D1 == 0 || obs.D2 == 0) continue;

            const NavData* nav = nullptr; double min_dt=1e9;
            for(const auto& n:navList) {
                if(n.sat==obs.sat && n.sys==obs.sys) {
                    double dt=timediff(obs.time, n.toe);
                    if(fabs(dt)<fabs(min_dt)) { min_dt=dt; nav=&n; }
                }
            }
            if(!nav) continue;

            double f1, f2;
            if(!get_sys_freqs(obs.sys, nav->freq_num, f1, f2)) continue;

            double lam1 = CLIGHT_GNSS / f1;
            double lam2 = CLIGHT_GNSS / f2;
            double rate1 = -obs.D1 * lam1;
            double rate2 = -obs.D2 * lam2;
            double obs_range_rate = calc_IF(rate1, rate2, f1, f2);

            double tau = obs.P1 / CLIGHT_GNSS;
            if(tau==0) tau = 0.075;
            GPSTime t_tx = obs.time; t_tx.sec -= tau;

            double rs[3]={0}, vs[3]={0}, dts=0.0, dts_drift=0.0;
            satPosVel(t_tx, *nav, rs, vs, &dts, &dts_drift);

            // =======================================================
            // [FIX] 速度解算也需要 TGD 吗？
            // 理论上 TGD 是常数，对多普勒(钟漂)无影响。
            // 但为了保持钟差 dts 的一致性，建议加上 (虽然 dts_drift 是导数，TGD 导数为 0)
            // 这里只修正 dts (用于 range rate model 中的 dts drift 其实不需要 TGD)
            // =======================================================

            double theta = OMEGA_E * tau;
            double cosT = cos(theta), sinT = sin(theta);
            double rx = rs[0]*cosT + rs[1]*sinT;
            double ry = -rs[0]*sinT + rs[1]*cosT;
            rs[0]=rx; rs[1]=ry;
            double vx = vs[0]*cosT + vs[1]*sinT;
            double vy = -vs[0]*sinT + vs[1]*cosT;
            vs[0]=vx; vs[1]=vy;

            double r = sqrt(pow(rs[0]-rr[0],2)+pow(rs[1]-rr[1],2)+pow(rs[2]-rr[2],2));
            double e[3] = {(rs[0]-rr[0])/r, (rs[1]-rr[1])/r, (rs[2]-rr[2])/r};

            double azel[2]={0}; sat_azel(rr, rs, azel);
            if(azel[1] < EL_MASK) continue;

            double rate_sat = e[0]*vs[0] + e[1]*vs[1] + e[2]*vs[2];
            double modeled_range_rate = rate_sat + CLIGHT_GNSS * dts_drift;

            double rate_rx = e[0]*vr[0] + e[1]*vr[1] + e[2]*vr[2];
            double rx_term = rate_rx + clk_drift[obs.sys];
            double b_val = obs_range_rate - (modeled_range_rate - rx_term);

            H.push_back(-e[0]); H.push_back(-e[1]); H.push_back(-e[2]);
            for(char s:syss) H.push_back(s==obs.sys?1.0:0.0);

            v.push_back(b_val);

            double w = 1.0;
            if(azel[1] > 0) w = sin(azel[1]) * sin(azel[1]);
            if (obs.sys == 'C' && obs.sat <= 5) w *= 0.1;
            W.push_back(w);
            nSat++;
        }

        if(nSat < m) return 0;

        vector<double> Q(m*m,0), b(m,0);
        for(int i=0;i<nSat;i++) {
            for(int j=0;j<m;j++) {
                for(int k=0;k<m;k++) Q[j*m+k]+=H[i*m+j]*W[i]*H[i*m+k];
                b[j]+=H[i*m+j]*W[i]*v[i];
            }
        }
        if(matinv(Q.data(), m)==-1) return 0;
        vector<double> dx(m);
        matmul("NN", m, m, 1, 1.0, Q.data(), b.data(), 0.0, dx.data());

        vr[0] += dx[0]; vr[1] += dx[1]; vr[2] += dx[2];
        for(int k=0;k<(int)syss.size();k++) clk_drift[syss[k]] += dx[3+k];

        if(sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]) < 1e-4) return nSat;
    }
    return nSat;
}

void spp_gnss_process(const std::vector<std::vector<ObsData>>& obsList,
                      const std::vector<NavData>& navList,
                      const double* refPos)
{
    cout << "==========================================================" << endl;
    cout << "   MYGNSS10 Multi-GNSS Dual-Freq IF SPP (FIXED TGD)" << endl;
    cout << "   Active Systems: ";
    if(USE_GPS) cout << "GPS ";
    if(USE_BDS) cout << "BDS ";
    if(USE_GAL) cout << "GAL ";
    if(USE_GLO) cout << "GLO ";
    cout << endl;
    cout << "==========================================================" << endl;
    cout << fixed << setprecision(3);

    cout << "GPST Time      Ns   PDOP    dE(m)    dN(m)    dU(m)      VE(m/s)  VN(m/s)  VU(m/s)" << endl;

    for (const auto& epochObs : obsList) {
        if (epochObs.empty()) continue;
        GPSTime t = epochObs[0].time;
        double ep[6]; time2epoch(t, ep);

        double current_hour = ep[3] + ep[4]/60.0 + ep[5]/3600.0;
        if (current_hour < 2.0 - 1e-9 || current_hour > 6.0 + 1e-9) continue;

        double rr[3] = {0};
        if(refPos[0]!=0) { rr[0]=refPos[0]; rr[1]=refPos[1]; rr[2]=refPos[2]; }

        map<char, double> clk;
        double pdop = 0.0;

        int ns = estPosition(epochObs, navList, rr, clk, &pdop);

        if(ns >= 4) {
            double vr[3] = {0};
            map<char, double> clk_drift;

            int nv = estVelocity(epochObs, navList, rr, vr, clk_drift);

            if(ep[5]>=59.95) {
                ep[5]=0; ep[4]++;
                if(ep[4]>=60){ep[4]=0; ep[3]++;}
            }
            char time_str[32];
            sprintf(time_str, "%02.0f:%02.0f:%04.1f", ep[3], ep[4], ep[5]);

            double enu_pos[3] = {0};
            if (refPos[0] != 0) xyz2enu(rr, refPos, enu_pos);
            else xyz2enu(rr, rr, enu_pos);

            double enu_vel[3] = {0};
            if (refPos[0] != 0) vel2enu(vr, refPos, enu_vel);
            else vel2enu(vr, rr, enu_vel);

            if(nv >= 4) {
                cout << time_str << "   "
                     << setw(2) << ns << "   "
                     << setw(5) << pdop << "   "
                     << setw(8) << enu_pos[0] << " "
                     << setw(8) << enu_pos[1] << " "
                     << setw(8) << enu_pos[2] << "      "
                     << setw(7) << enu_vel[0] << "  "
                     << setw(7) << enu_vel[1] << "  "
                     << setw(7) << enu_vel[2] << endl;
            }
        }
    }
}