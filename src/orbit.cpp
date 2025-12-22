/*
 * src/orbit.cpp - 全系统卫星轨道计算 (GPS/BDS/GAL/GLO)
 * 核心功能:
 * 1. 开普勒轨道 (GPS/BDS/GAL)
 * 2. 数值积分轨道 (GLONASS)
 * 3. 卫星速度与钟漂计算 (支持多普勒测速)
 */

#include "orbit.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

// --- 物理常数定义 ---
const double PI = 3.1415926535897932;

// 地球引力常数 (GM)
const double GM_GPS = 3.9860050E14;
const double GM_BDS = 3.986004418E14;
const double GM_GAL = 3.986004418E14;
const double GM_GLO = 3.986004418E14;

// 地球自转角速度 (Omega_e)
const double OMEGA_GPS = 7.2921151467E-5;
const double OMEGA_BDS = 7.2921150E-5;
const double OMEGA_GAL = 7.2921151467E-5;
const double OMEGA_GLO = 7.2921150E-5;

// GLONASS 特定常数
const double RE_GLO = 6378136.0;   // 地球长半轴
const double J2_GLO = 1.0826257E-3;// 二阶带谐系数

// BDS GEO 旋转角 (-5度)
const double BDS_GEO_ANGLE = -5.0 * PI / 180.0;

// --- 辅助函数 ---

// 判断是否为 BDS GEO (C01-C05)
static bool is_bds_geo(const NavData& nav) {
    if (nav.sys != 'C') return false;
    int prn = nav.sat;
    if (prn > 200) prn -= 200; // 处理编号偏移
    return (prn >= 1 && prn <= 5);
}

// --------------------------------------------------------------------------
// 模块 1: 开普勒轨道计算核心 (GPS / BDS / Galileo)
// --------------------------------------------------------------------------
// 内部函数：仅计算位置和钟差
static void satPosKeplerCore(GPSTime tTx, const NavData& nav, double* rs, double* dts) {
    // 1. 选择常数
    double GM_VAL, OMEGA_VAL;
    if (nav.sys == 'C') { GM_VAL = GM_BDS; OMEGA_VAL = OMEGA_BDS; }
    else if (nav.sys == 'E') { GM_VAL = GM_GAL; OMEGA_VAL = OMEGA_GAL; }
    else { GM_VAL = GM_GPS; OMEGA_VAL = OMEGA_GPS; }

    // 2. 时间归化
    double tk = timediff(tTx, nav.toe);

    // 3. 基础参数
    double n0 = sqrt(GM_VAL) / pow(nav.sqrtA, 3);
    double n = n0 + nav.delta_n;
    double M = nav.M0 + n * tk;

    // 4. 开普勒方程迭代
    double E = M, Eold;
    for (int i = 0; i < 30; i++) {
        Eold = E;
        E = M + nav.e * sin(E);
        if (fabs(E - Eold) < 1e-13) break;
    }

    // 5. 真近点角 & 升交角距
    double sinE = sin(E), cosE = cos(E);
    double nu = atan2(sqrt(1.0 - nav.e*nav.e)*sinE, cosE - nav.e);
    double phi = nu + nav.omega;

    // 6. 摄动改正
    double sin2phi = sin(2.0 * phi), cos2phi = cos(2.0 * phi);
    double u = phi + nav.cus * sin2phi + nav.cuc * cos2phi;
    double r = nav.sqrtA*nav.sqrtA * (1.0 - nav.e * cosE) + nav.crs * sin2phi + nav.crc * cos2phi;
    double i = nav.i0 + nav.IDot * tk + nav.cis * sin2phi + nav.cic * cos2phi;

    // 7. 轨道平面坐标
    double x_prime = r * cos(u);
    double y_prime = r * sin(u);

    // 8. 转换至 ECEF (含 BDS GEO 特殊处理)
    double Omega = nav.Omega0 + (nav.OmegaDot - OMEGA_VAL) * tk - OMEGA_VAL * nav.toe.sec;
    double sinO = sin(Omega), cosO = cos(Omega);
    double sini = sin(i), cosi = cos(i);

    // 惯性系/中间坐标
    double x_tmp = x_prime * cosO - y_prime * cosi * sinO;
    double y_tmp = x_prime * sinO + y_prime * cosi * cosO;
    double z_tmp = y_prime * sini;

    if (is_bds_geo(nav)) {
        // BDS GEO: R_x(-5 deg)
        double ang = BDS_GEO_ANGLE;
        double sinA = sin(ang), cosA = cos(ang);
        rs[0] = x_tmp;
        rs[1] = y_tmp * cosA - z_tmp * sinA;
        rs[2] = y_tmp * sinA + z_tmp * cosA;
    } else {
        rs[0] = x_tmp;
        rs[1] = y_tmp;
        rs[2] = z_tmp;
    }

    // 9. 钟差 (含相对论)
    if (dts) {
        *dts = nav.a0 + nav.a1 * tk + nav.a2 * tk * tk;
        *dts -= 2.0 * sqrt(GM_VAL) * nav.sqrtA * nav.e * sinE / (CLIGHT * CLIGHT);
    }
}

// --------------------------------------------------------------------------
// 模块 2: GLONASS 数值积分 (Runge-Kutta 4)
// --------------------------------------------------------------------------

// GLONASS 微分方程: dy/dt = f(t, y)
// y = [x, y, z, vx, vy, vz]
static void glonass_deq(const double *y, double *yp, const double *acc) {
    double r2 = y[0]*y[0] + y[1]*y[1] + y[2]*y[2];
    double r3 = r2 * sqrt(r2);
    double r5 = r2 * r2 * sqrt(r2);

    // J2 项与引力项
    double a = 1.5 * J2_GLO * GM_GLO * RE_GLO * RE_GLO / r5;
    double b = 5.0 * y[2] * y[2] / r2;
    double c = -GM_GLO / r3 + a * (1.0 - b);

    // Pos 导数 = Vel
    yp[0] = y[3]; yp[1] = y[4]; yp[2] = y[5];

    // Vel 导数 = Acc
    yp[3] = c * y[0] + OMEGA_GLO * OMEGA_GLO * y[0] + 2.0 * OMEGA_GLO * y[4] + acc[0];
    yp[4] = c * y[1] + OMEGA_GLO * OMEGA_GLO * y[1] - 2.0 * OMEGA_GLO * y[3] + acc[1];
    yp[5] = (c - 2.0 * a) * y[2] + acc[2];
}

// RK4 积分器 (支持输出速度)
static void glonass_orbit(double t, const double *x, double *rs, double *vs) {
    // x[0-2]: pos, x[3-5]: vel, x[6-8]: acc
    double y[6], k1[6], k2[6], k3[6], k4[6], w[6];
    double acc[3] = {x[6], x[7], x[8]};

    for(int i=0; i<6; i++) y[i] = x[i];

    double h = (t < 0.0) ? -30.0 : 30.0;

    for (double tt = 0.0; fabs(tt) < fabs(t) + 1e-9; tt += h) {
        if (fabs(t - tt) < fabs(h)) h = t - tt;
        if (fabs(h) < 1e-9) break;

        glonass_deq(y, k1, acc);
        for(int i=0; i<6; i++) w[i] = y[i] + k1[i] * h / 2.0;
        glonass_deq(w, k2, acc);
        for(int i=0; i<6; i++) w[i] = y[i] + k2[i] * h / 2.0;
        glonass_deq(w, k3, acc);
        for(int i=0; i<6; i++) w[i] = y[i] + k3[i] * h;
        glonass_deq(w, k4, acc);
        for(int i=0; i<6; i++) y[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) * h / 6.0;
    }

    if (rs) { rs[0] = y[0]; rs[1] = y[1]; rs[2] = y[2]; }
    if (vs) { vs[0] = y[3]; vs[1] = y[4]; vs[2] = y[5]; }
}

// --------------------------------------------------------------------------
// 统一入口
// --------------------------------------------------------------------------

// 接口1: 仅定位 (保留原有功能)
void satPos(GPSTime tTx, const NavData& nav, double* rs, double* dts) {
    if (nav.sys == 'R') {
        double x[9];
        x[0] = nav.M0;     x[1] = nav.e;      x[2] = nav.sqrtA;  // Pos
        x[3] = nav.Omega0; x[4] = nav.i0;     x[5] = nav.omega;  // Vel
        x[6] = nav.cuc;    x[7] = nav.cus;    x[8] = nav.crc;    // Acc

        double t = timediff(tTx, nav.toe);
        glonass_orbit(t, x, rs, nullptr);
        if (dts) *dts = nav.a0 + nav.a1 * t;
    } else {
        satPosKeplerCore(tTx, nav, rs, dts);
    }
}

// 接口2: 定位 + 测速 (新增功能)
void satPosVel(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift) {
    if (nav.sys == 'R') {
        // --- GLONASS 处理 ---
        double x[9];
        x[0] = nav.M0;     x[1] = nav.e;      x[2] = nav.sqrtA;
        x[3] = nav.Omega0; x[4] = nav.i0;     x[5] = nav.omega;
        x[6] = nav.cuc;    x[7] = nav.cus;    x[8] = nav.crc;

        double t = timediff(tTx, nav.toe);
        // 直接从积分结果获取位置和速度
        glonass_orbit(t, x, rs, vs);

        // 钟差与钟漂
        if (dts) *dts = nav.a0 + nav.a1 * t;
        if (dts_drift) *dts_drift = nav.a1; // 线性钟漂即 gamma_n
    } else {
        // --- Kepler 处理 (GPS/BDS/GAL) ---
        // 使用中心差分法求速度 (Robust way)
        double dt = 0.001; // 1ms 步长
        double rs1[3], rs2[3], dts1, dts2;

        GPSTime t1 = tTx; t1.sec -= dt;
        GPSTime t2 = tTx; t2.sec += dt;

        // 计算中心时刻位置钟差
        satPosKeplerCore(tTx, nav, rs, dts);

        // 计算前后时刻位置钟差
        satPosKeplerCore(t1, nav, rs1, &dts1);
        satPosKeplerCore(t2, nav, rs2, &dts2);

        // 差分求速度
        if (vs) {
            vs[0] = (rs2[0] - rs1[0]) / (2.0 * dt);
            vs[1] = (rs2[1] - rs1[1]) / (2.0 * dt);
            vs[2] = (rs2[2] - rs1[2]) / (2.0 * dt);
        }

        // 差分求钟漂 (自动包含了相对论效应的变化率)
        if (dts_drift) {
            *dts_drift = (dts2 - dts1) / (2.0 * dt);
        }
    }
}