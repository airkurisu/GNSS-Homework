/*
 * src/spp.cpp - GPS单系统 SPP定位与测速实现
 */

#include "spp.h"
#include "gtime.h"
#include "matrix.h"
#include "orbit.h"
#include "tropo.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdio>

using namespace std;

// --- 常量定义 ---
//const double CLIGHT = 299792458.0; // 在 orbit.h 或其他地方可能已定义，注意避免重定义
const double PI = 3.1415926535897932;
const double OMEGA_E = 7.2921151467E-5;   // 地球自转角速度
const double EL_MASK = 10.0 * PI / 180.0; // 截止高度角 10度

// 椭球参数 (WGS84)
const double RE_WGS84 = 6378137.0;
const double FE_WGS84 = 1.0 / 298.257223563;

// 频率定义 (GPS)
const double FREQ_L1 = 1.57542E9;
const double FREQ_L2 = 1.22760E9;
const double LAMBDA_L1 = CLIGHT / FREQ_L1;
const double LAMBDA_L2 = CLIGHT / FREQ_L2;

// 无电离层组合系数
const double ALPHA = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
const double BETA  = (FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);

// --- 辅助工具函数 ---

// ECEF (XYZ) 转 Geodetic (Lat, Lon, H)
// pos: {lat, lon, h} (rad, rad, m)
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
    double sinP = sin(lat), cosP = cos(lat);
    double sinL = sin(lon), cosL = cos(lon);
    double dx = r[0] - ref[0];
    double dy = r[1] - ref[1];
    double dz = r[2] - ref[2];
    enu[0] = -sinL * dx + cosL * dy;
    enu[1] = -sinP * cosL * dx - sinP * sinL * dy + cosP * dz;
    enu[2] =  cosP * cosL * dx + cosP * sinL * dy + sinP * dz;
}

static void vel2enu(const double *v_xyz, const double *ref_pos, double *v_enu) {
    double r_dummy[3] = { ref_pos[0] + v_xyz[0], ref_pos[1] + v_xyz[1], ref_pos[2] + v_xyz[2] };
    xyz2enu(r_dummy, ref_pos, v_enu);
}

static void sat_azel(const double *r, const double *rs, double *azel) {
    double enu[3];
    xyz2enu(rs, r, enu);
    double dist = sqrt(enu[0]*enu[0] + enu[1]*enu[1] + enu[2]*enu[2]);
    if (dist < 1e-3) { azel[0]=0; azel[1]=0; return; }
    azel[1] = asin(enu[2] / dist);
    azel[0] = atan2(enu[0], enu[1]);
    if (azel[0] < 0) azel[0] += 2*PI;
}

// ------------------------------------------------------------------------------------------------
// 核心模块 1: 伪距单点定位
// ------------------------------------------------------------------------------------------------
static int estPosition(const vector<ObsData>& obsList, const vector<NavData>& navList,
                       double* x, double* dtr, double* pdop) {
    int nSat = 0;
    for (int iter = 0; iter < 10; iter++) {
        vector<double> H, v, W;
        nSat = 0;

        // [Fix] 在调用对流层模型前，计算当前坐标的大地坐标 (BLH)
        double pos[3] = {0};
        ecef2pos(x, pos);

        for (const auto& obs : obsList) {
            if (obs.sys != 'G') continue;
            if (obs.P1 == 0.0 || obs.P2 == 0.0) continue;

            const NavData* nav = nullptr;
            double min_dt = 1e9;
            for (const auto& n : navList) {
                if (n.sat == obs.sat && n.sys == obs.sys) {
                    double dt = timediff(obs.time, n.toe);
                    if (fabs(dt) < fabs(min_dt)) { min_dt = dt; nav = &n; }
                }
            }
            if (!nav || fabs(min_dt) > 7200.0) continue;

            GPSTime t_tx = obs.time;
            t_tx.sec -= obs.P1 / CLIGHT;

            double rs[3], dts;
            satPos(t_tx, *nav, rs, &dts);

            double r_dist = sqrt(pow(rs[0]-x[0],2) + pow(rs[1]-x[1],2) + pow(rs[2]-x[2],2));
            double tau = r_dist / CLIGHT;
            double theta = OMEGA_E * tau;
            double rx = rs[0]*cos(theta) + rs[1]*sin(theta);
            double ry = -rs[0]*sin(theta) + rs[1]*cos(theta);
            rs[0] = rx; rs[1] = ry;

            double r = sqrt(pow(rs[0]-x[0],2) + pow(rs[1]-x[1],2) + pow(rs[2]-x[2],2));
            double ex = (x[0] - rs[0]) / r;
            double ey = (x[1] - rs[1]) / r;
            double ez = (x[2] - rs[2]) / r;

            double azel[2] = {0, 0};
            if (iter == 0 && sqrt(x[0]*x[0]+x[1]*x[1]) < 1.0) azel[1] = PI/2.0;
            else sat_azel(x, rs, azel);

            if (azel[1] < EL_MASK) continue;

            double trop = 0.0, trop_var = 0.0;
            // [Fix] 传入 BLH 坐标 pos，而不是 ECEF 坐标 x
            if (iter > 0) trop_model_prec(obs.time, pos, azel, &trop, &trop_var);

            double P_IF = ALPHA * obs.P1 - BETA * obs.P2;
            double P_comp = r + *dtr - CLIGHT * dts + trop;

            H.push_back(ex); H.push_back(ey); H.push_back(ez); H.push_back(1.0);
            v.push_back(P_IF - P_comp);
            W.push_back(sin(azel[1])*sin(azel[1]));
            nSat++;
        }

        if (nSat < 4) return 0;

        int m = 4;
        vector<double> Q(m*m, 0.0), Q_geom(m*m, 0.0), b(m, 0.0);

        for(int i=0; i<nSat; i++) {
            double w = W[i];
            for(int j=0; j<4; j++) {
                for(int k=0; k<4; k++) {
                    Q[j*4+k] += H[i*4+j] * w * H[i*4+k];
                    Q_geom[j*4+k] += H[i*4+j] * 1.0 * H[i*4+k];
                }
                b[j] += H[i*4+j] * w * v[i];
            }
        }

        if (matinv(Q.data(), m) == -1) return 0;
        double dx[4] = {0};
        matmul("NN", 4, 4, 1, 1.0, Q.data(), b.data(), 0.0, dx);

        x[0] += dx[0]; x[1] += dx[1]; x[2] += dx[2]; *dtr += dx[3];

        if (matinv(Q_geom.data(), m) != -1) *pdop = sqrt(Q_geom[0] + Q_geom[5] + Q_geom[10]);
        else *pdop = 0.0;

        if (sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) < 1e-4) return nSat;
    }
    return 0;
}

// ------------------------------------------------------------------------------------------------
// 核心模块 2: 多普勒测速
// ------------------------------------------------------------------------------------------------
static int estVelocity(const vector<ObsData>& obsList, const vector<NavData>& navList,
                       const double* rr, double* vr, double* dtr_drift) {
    int nSat = 0;
    for (int iter = 0; iter < 5; iter++) {
        vector<double> H, v, W;
        nSat = 0;
        for (const auto& obs : obsList) {
            if (obs.sys != 'G') continue;
            if (obs.D1 == 0.0 || obs.D2 == 0.0) continue;

            const NavData* nav = nullptr;
            double min_dt = 1e9;
            for (const auto& n : navList) {
                if (n.sat == obs.sat && n.sys == obs.sys) {
                    double dt = timediff(obs.time, n.toe);
                    if (fabs(dt) < fabs(min_dt)) { min_dt = dt; nav = &n; }
                }
            }
            if (!nav) continue;

            GPSTime t_tx = obs.time; t_tx.sec -= 0.075;
            double rs[3], vs[3], dts, dts_drift;
            satPosVel(t_tx, *nav, rs, vs, &dts, &dts_drift);

            double azel[2] = {0, 0};
            sat_azel(rr, rs, azel);
            if (azel[1] < EL_MASK) continue;

            double r = sqrt(pow(rs[0]-rr[0],2) + pow(rs[1]-rr[1],2) + pow(rs[2]-rr[2],2));
            double e[3] = { (rs[0]-rr[0])/r, (rs[1]-rr[1])/r, (rs[2]-rr[2])/r };

            double rate_L1 = -obs.D1 * LAMBDA_L1;
            double rate_L2 = -obs.D2 * LAMBDA_L2;
            double obs_rate_IF = ALPHA * rate_L1 - BETA * rate_L2;
            double rate_sat_proj = e[0]*vs[0] + e[1]*vs[1] + e[2]*vs[2];
            double rate_pred = rate_sat_proj - (e[0]*vr[0] + e[1]*vr[1] + e[2]*vr[2])
                               - CLIGHT * dts_drift + *dtr_drift;

            v.push_back(obs_rate_IF - rate_pred);
            H.push_back(-e[0]); H.push_back(-e[1]); H.push_back(-e[2]); H.push_back(1.0);
            W.push_back(1.0);
            nSat++;
        }

        if (nSat < 4) return 0;
        int m = 4;
        vector<double> Q(m*m, 0.0), b(m, 0.0);
        for(int i=0; i<nSat; i++) {
            for(int j=0; j<4; j++) {
                for(int k=0; k<4; k++) Q[j*4+k] += H[i*4+j] * W[i] * H[i*4+k];
                b[j] += H[i*4+j] * W[i] * v[i];
            }
        }
        if (matinv(Q.data(), m) == -1) return 0;
        double dx[4] = {0};
        matmul("NN", 4, 4, 1, 1.0, Q.data(), b.data(), 0.0, dx);

        vr[0] += dx[0]; vr[1] += dx[1]; vr[2] += dx[2]; *dtr_drift += dx[3];
        if (sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) < 1e-4) return nSat;
    }
    return nSat;
}

// ------------------------------------------------------------------------------------------------
// 主流程控制
// ------------------------------------------------------------------------------------------------
void spp_process(const std::vector<std::vector<ObsData>>& obsList,
                 const std::vector<NavData>& navList,
                 const double* refPos)
{
    cout << "==========================================================" << endl;
    cout << "   MYGNSS10 SPP Processing (GPS Only)" << endl;
    cout << "   Time Filter: 02:00 - 06:00 (GPST)" << endl;
    cout << "   PDOP: Standard Geometric (Unweighted)" << endl;
    cout << "==========================================================" << endl;
    cout << fixed << setprecision(3);
    // [Update] 输出表头增加 ZTD(m)
    cout << "GPST Time    Ns   PDOP       dE(m)    dN(m)    dU(m)   VE(m/s)  VN(m/s)  VU(m/s)" << endl;

    for (const auto& epochObs : obsList) {
        if (epochObs.empty()) continue;
        GPSTime t_obs = epochObs[0].time;

        double ep_check[6];
        time2epoch(t_obs, ep_check);
        double current_hour = ep_check[3] + ep_check[4]/60.0 + ep_check[5]/3600.0;
        if (current_hour < 2.0 - 1e-9 || current_hour > 6.0 + 1e-9) continue;

        double rr[3] = {0};
        double dtr = 0.0;
        double pdop = 0.0;
        if (refPos[0] != 0) { rr[0]=refPos[0]; rr[1]=refPos[1]; rr[2]=refPos[2]; }

        int nSatPos = estPosition(epochObs, navList, rr, &dtr, &pdop);

        if (nSatPos >= 4) {
            double vr[3] = {0};
            double dtr_drift = 0.0;
            estVelocity(epochObs, navList, rr, vr, &dtr_drift);

            // --- [Update] 计算当前历元的天顶对流层延迟 (ZTD) ---
            double ztd = 0.0, dummy_var = 0.0;
            double pos_blh[3] = {0};
            double azel_zenith[2] = {0.0, PI / 2.0}; // 方位角0，高度角90度

            ecef2pos(rr, pos_blh); // 转为 BLH
            // 计算 ZTD
            trop_model_prec(t_obs, pos_blh, azel_zenith, &ztd, &dummy_var);

            // --- 修复时间显示逻辑 ---
            double ep[6];
            time2epoch(t_obs, ep);
            if (ep[5] >= 59.95) {
                ep[5] = 0.0;
                ep[4] += 1.0;
                if (ep[4] >= 60.0) { ep[4] = 0.0; ep[3] += 1.0; }
            }

            char time_str[32];
            sprintf(time_str, "%02.0f:%02.0f:%04.1f", ep[3], ep[4], ep[5]);

            double enu_pos[3] = {0};
            if (refPos[0] != 0) xyz2enu(rr, refPos, enu_pos);
            else xyz2enu(rr, rr, enu_pos);

            double enu_vel[3] = {0};
            vel2enu(vr, rr, enu_vel);

            // [Update] 输出 ZTD
            cout << time_str << "   "
            << setw(2) << nSatPos << "   "
            << setw(5) << pdop << "   "
            << setw(8) << enu_pos[0] << " "
            << setw(8) << enu_pos[1] << " "
            << setw(8) << enu_pos[2] << " "
            << setw(7) << enu_vel[0] << "  "
            << setw(7) << enu_vel[1] << "  "
            << setw(7) << enu_vel[2] << endl;
        }
    }
}