#include "orbit.h"
#include <cmath>
#include <iostream>
#include <cstring>

using namespace std;

// --- 物理常数 ---
const double CLIGHT_VAL = 299792458.0;
const double PI = 3.1415926535897932;

const double GM_GPS = 3.9860050E14;
const double GM_BDS = 3.986004418E14;
const double GM_GAL = 3.986004418E14;
const double GM_GLO = 3.9860044E14;

const double OMEGA_GPS = 7.2921151467E-5;
const double OMEGA_BDS = 7.2921150E-5;
const double OMEGA_GAL = 7.2921151467E-5;
const double OMEGA_GLO = 7.2921150E-5;

// GLONASS 特定参数
const double RE_GLO = 6378136.0;
const double J2_GLO = 1.0826257E-3;

// BDS GEO 5度倾角参数
const double SIN_5 = -0.0871557427476582;
const double COS_5 =  0.9961946980917456;

// --------------------------------------------------------------------------
// 模块 1: 开普勒轨道 + 解析速度 (GPS / BDS / Galileo)
// --------------------------------------------------------------------------
static void satPosKepler(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift) {
    double GM_VAL, OMEGA_VAL;
    if (nav.sys == 'C') { GM_VAL = GM_BDS; OMEGA_VAL = OMEGA_BDS; }
    else if (nav.sys == 'E') { GM_VAL = GM_GAL; OMEGA_VAL = OMEGA_GAL; }
    else { GM_VAL = GM_GPS; OMEGA_VAL = OMEGA_GPS; }

    // 1. 计算归化时间 tk
    double tk = timediff(tTx, nav.toe);

    // [BDS 时间修正] BDT -> GPST (14s)
    // BDT was already converted to GPST by the reader.

    // 2. 平近点角 M
    double A = nav.sqrtA * nav.sqrtA;
    double n0 = sqrt(GM_VAL / (A * A * A));
    double n = n0 + nav.delta_n;
    double M = nav.M0 + n * tk;

    // 3. 偏近点角 E (迭代求解)
    double E = M, Eold;
    for (int i = 0; i < 30; i++) {
        Eold = E;
        E = M + nav.e * sin(E);
        if (fabs(E - Eold) < 1e-13) break;
    }

    double sinE = sin(E);
    double cosE = cos(E);
    double Edot = n / (1.0 - nav.e * cosE); 

    // 4. 真近点角 nu
    double phi = atan2(sqrt(1.0 - nav.e * nav.e) * sinE, cosE - nav.e) + nav.omega;
    double sin2phi = sin(2.0 * phi);
    double cos2phi = cos(2.0 * phi);

    // 5. 摄动改正
    double du = nav.cus * sin2phi + nav.cuc * cos2phi;
    double dr = nav.crs * sin2phi + nav.crc * cos2phi;
    double di = nav.cis * sin2phi + nav.cic * cos2phi;

    double u = phi + du;
    double r = A * (1.0 - nav.e * cosE) + dr;
    double i = nav.i0 + nav.IDot * tk + di;

    // 6. 轨道平面坐标
    double x_prime = r * cos(u);
    double y_prime = r * sin(u);

    // 导数计算
    double nudot = Edot * sqrt(1.0 - nav.e * nav.e) / (1.0 - nav.e * cosE);
    double udot = nudot + 2.0 * nudot * (nav.cus * cos2phi - nav.cuc * sin2phi);
    double rdot = A * nav.e * sinE * Edot + 2.0 * nudot * (nav.crs * cos2phi - nav.crc * sin2phi);
    double idot = nav.IDot + 2.0 * nudot * (nav.cis * cos2phi - nav.cic * sin2phi);
    
    double vx_prime = rdot * cos(u) - r * sin(u) * udot;
    double vy_prime = rdot * sin(u) + r * cos(u) * udot;

    // 7. 坐标转换 (区分 BDS GEO)
    
    bool is_bds_geo = (nav.sys == 'C' && (nav.sat <= 5 || nav.sat >= 59));
    
    double x_final, y_final, z_final;
    double vx_final, vy_final, vz_final;
    
    // [Fix] 提前声明变量
    double Omega, dOmega;

    if (is_bds_geo) {
        // --- BDS GEO 特殊处理 (ICD 5.1.4) ---
        
        // A. 计算升交点经度 Omega (不减 Omega_e * tk)
        Omega = nav.Omega0 + nav.OmegaDot * tk - OMEGA_VAL * nav.toeSow;
        dOmega = nav.OmegaDot; 

        double sinO = sin(Omega), cosO = cos(Omega);
        double sini = sin(i), cosi = cos(i);

        // B. 计算在 "GEO 惯性系" 中的坐标
        double xg = x_prime * cosO - y_prime * cosi * sinO;
        double yg = x_prime * sinO + y_prime * cosi * cosO;
        double zg = y_prime * sini;

        // C. 计算 "GEO 惯性系" 中的速度
        double vxg = vx_prime * cosO - vy_prime * cosi * sinO 
                   - (x_prime * sinO + y_prime * cosi * cosO) * dOmega 
                   + y_prime * sinO * sini * idot;

        double vyg = vx_prime * sinO + vy_prime * cosi * cosO 
                   + (x_prime * cosO - y_prime * cosi * sinO) * dOmega 
                   - y_prime * cosO * sini * idot;

        double vzg = vy_prime * sini + y_prime * cosi * idot;

        // D. 坐标转换到 ECEF
        // Rx(-5)
        double rx = xg;
        double ry = yg * COS_5 + zg * SIN_5;
        double rz = -yg * SIN_5 + zg * COS_5;

        // Rz(omega*tk)
        double alpha = OMEGA_VAL * tk;
        double sinA = sin(alpha);
        double cosA = cos(alpha);

        x_final = rx * cosA + ry * sinA;
        y_final = -rx * sinA + ry * cosA;
        z_final = rz;

        // E. 速度转换
        double vrx = vxg;
        double vry = vyg * COS_5 + vzg * SIN_5;
        double vrz = -vyg * SIN_5 + vzg * COS_5;

        double dPx = rx * (-sinA * OMEGA_VAL) + ry * (cosA * OMEGA_VAL);
        double dPy = -rx * (cosA * OMEGA_VAL) + ry * (-sinA * OMEGA_VAL);
        double dPz = 0.0;

        double dVx = vrx * cosA + vry * sinA;
        double dVy = -vrx * sinA + vry * cosA;
        double dVz = vrz;

        vx_final = dPx + dVx;
        vy_final = dPy + dVy;
        vz_final = dPz + dVz;

    } else {
        // --- MEO / IGSO (标准公式) ---
        Omega = nav.Omega0 + (nav.OmegaDot - OMEGA_VAL) * tk - OMEGA_VAL * nav.toeSow;
        dOmega = nav.OmegaDot - OMEGA_VAL;

        double sinO = sin(Omega), cosO = cos(Omega);
        double sini = sin(i), cosi = cos(i);

        x_final = x_prime * cosO - y_prime * cosi * sinO;
        y_final = x_prime * sinO + y_prime * cosi * cosO;
        z_final = y_prime * sini;

        vx_final = (vx_prime * cosO - vy_prime * cosi * sinO) -
                   (x_prime * sinO + y_prime * cosi * cosO) * dOmega +
                   y_prime * sinO * sini * idot;

        vy_final = (vx_prime * sinO + vy_prime * cosi * cosO) +
                   (x_prime * cosO - y_prime * cosi * sinO) * dOmega -
                   y_prime * cosO * sini * idot;

        vz_final = vy_prime * sini + y_prime * cosi * idot;
    }

    // 8. 赋值
    rs[0] = x_final; rs[1] = y_final; rs[2] = z_final;
    if (vs) { vs[0] = vx_final; vs[1] = vy_final; vs[2] = vz_final; }

    const double tc = timediff(tTx, nav.toc);

    // 9. 卫星钟差
    if (dts) {
        *dts = nav.a0 + nav.a1 * tc + nav.a2 * tc * tc;
        *dts -= 2.0 * sqrt(GM_VAL * A) * nav.e * sinE / (CLIGHT_VAL * CLIGHT_VAL);
    }

    // 10. 卫星钟漂
    if (dts_drift) {
        *dts_drift = nav.a1 + 2.0 * nav.a2 * tc
            - 2.0 * sqrt(GM_VAL * A) * nav.e * cosE * Edot / (CLIGHT_VAL * CLIGHT_VAL);
    }
}

// --------------------------------------------------------------------------
// 模块 2: GLONASS 轨道积分 (保持修正)
// --------------------------------------------------------------------------
static void glonass_deq(const double *x, double *xdot, const double *acc) {
    double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    double r3 = r2 * sqrt(r2);
    double r5 = r2 * r2 * sqrt(r2);

    double a = 1.5 * J2_GLO * GM_GLO * RE_GLO * RE_GLO / r5;
    double b = 5.0 * x[2] * x[2] / r2;
    double c = -GM_GLO / r3 - a * (1.0 - b);

    xdot[0] = x[3]; xdot[1] = x[4]; xdot[2] = x[5];

    double omg2 = OMEGA_GLO * OMEGA_GLO;
    xdot[3] = (c + omg2) * x[0] + 2.0 * OMEGA_GLO * x[4] + acc[0];
    xdot[4] = (c + omg2) * x[1] - 2.0 * OMEGA_GLO * x[3] + acc[1];
    xdot[5] = (c - 2.0 * a) * x[2] + acc[2];
}

static void glonass_orbit(double t, double *x, const double *acc) {
    double k1[6], k2[6], k3[6], k4[6], w[6];
    double step = (t < 0.0) ? -30.0 : 30.0;
    double dt = t;

    while (fabs(dt) > 1e-9) {
        if (fabs(dt) < fabs(step)) step = dt;

        glonass_deq(x, k1, acc); for(int i=0; i<6; i++) w[i] = x[i] + k1[i] * step / 2.0;
        glonass_deq(w, k2, acc); for(int i=0; i<6; i++) w[i] = x[i] + k2[i] * step / 2.0;
        glonass_deq(w, k3, acc); for(int i=0; i<6; i++) w[i] = x[i] + k3[i] * step;
        glonass_deq(w, k4, acc); for(int i=0; i<6; i++) x[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) * step / 6.0;

        dt -= step;
    }
}

void satPosVel(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift) {
    if (nav.sys == 'R') {
        double x[6] = { nav.gloPos[0], nav.gloPos[1], nav.gloPos[2], nav.gloVel[0], nav.gloVel[1], nav.gloVel[2] };
        double acc[3] = { nav.gloAcc[0], nav.gloAcc[1], nav.gloAcc[2] };
        double t = timediff(tTx, nav.toe);
        //t -= 18.0; // 闰秒修正
        glonass_orbit(t, x, acc);
        if (rs) { rs[0]=x[0]; rs[1]=x[1]; rs[2]=x[2]; }
        if (vs) { vs[0]=x[3]; vs[1]=x[4]; vs[2]=x[5]; }
        if (dts) *dts = nav.a0 + nav.a1 * t;
        if (dts_drift) *dts_drift = nav.a1;
    } else {
        satPosKepler(tTx, nav, rs, vs, dts, dts_drift);
    }
}

void satPos(GPSTime tTx, const NavData& nav, double* rs, double* dts) {
    satPosVel(tTx, nav, rs, nullptr, dts, nullptr);
}
