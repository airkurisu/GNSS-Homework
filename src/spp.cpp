/*
 * src/spp.cpp - SPP定位与测速实现
 * 更新：增加时间过滤 (2:00 - 6:00)
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
#include <cstdio> // 确保包含 printf/sprintf

using namespace std;

// --- 常量定义 ---
const double PI = 3.1415926535897932;
const double OMEGA_E = 7.2921151467E-5; // 地球自转角速度

// 频率定义 (GPS)
const double FREQ_L1 = 1.57542E9;
const double FREQ_L2 = 1.22760E9;
const double LAMBDA_L1 = CLIGHT / FREQ_L1; // L1波长 (~0.19m)
const double LAMBDA_L2 = CLIGHT / FREQ_L2; // L2波长 (~0.24m)

// 无电离层组合系数
const double ALPHA = (FREQ_L1 * FREQ_L1) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);
const double BETA  = (FREQ_L2 * FREQ_L2) / (FREQ_L1 * FREQ_L1 - FREQ_L2 * FREQ_L2);

// --- 辅助工具函数 ---

// 坐标转换 XYZ -> ENU
static void xyz2enu(const double *r, const double *ref, double *enu) {
    double lat = atan2(ref[2], sqrt(ref[0]*ref[0]+ref[1]*ref[1]));
    double lon = atan2(ref[1], ref[0]);
    double sinP = sin(lat), cosP = cos(lat);
    double sinL = sin(lon), cosL = cos(lon);

    double dx = r[0] - ref[0];
    double dy = r[1] - ref[1];
    double dz = r[2] - ref[2];

    enu[0] = -sinL * dx + cosL * dy;                      // E
    enu[1] = -sinP * cosL * dx - sinP * sinL * dy + cosP * dz; // N
    enu[2] =  cosP * cosL * dx + cosP * sinL * dy + sinP * dz; // U
}

// 速度转换 XYZ -> ENU (旋转矩阵与位置相同)
static void vel2enu(const double *v_xyz, const double *ref_pos, double *v_enu) {
    double r_dummy[3] = { ref_pos[0] + v_xyz[0], ref_pos[1] + v_xyz[1], ref_pos[2] + v_xyz[2] };
    xyz2enu(r_dummy, ref_pos, v_enu);
}

// 计算卫星方位角/高度角
static void sat_azel(const double *r, const double *rs, double *azel) {
    double enu[3];
    xyz2enu(rs, r, enu);
    double dist = sqrt(enu[0]*enu[0] + enu[1]*enu[1] + enu[2]*enu[2]);
    if (dist < 1e-3) { azel[0]=0; azel[1]=0; return; }
    azel[1] = asin(enu[2] / dist); // El (rad)
    azel[0] = atan2(enu[0], enu[1]); // Az (rad)
    if (azel[0] < 0) azel[0] += 2*PI;
}

// ------------------------------------------------------------------------------------------------
// 核心模块 1: 伪距单点定位 (Positioning)
// ------------------------------------------------------------------------------------------------
static int estPosition(const vector<ObsData>& obsList, const vector<NavData>& navList,
                       double* x, double* dtr, double* pdop) {
    int nSat = 0;

    // 迭代求解
    for (int iter = 0; iter < 10; iter++) {
        vector<double> H, v, W;
        nSat = 0;

        for (const auto& obs : obsList) {
            // 这里以 GPS 为例，若要扩展多系统需修改此处逻辑
            if (obs.sys != 'G') continue;
            if (obs.P1 == 0.0 || obs.P2 == 0.0) continue;

            // 1. 匹配星历
            const NavData* nav = nullptr;
            double min_dt = 1e9;
            for (const auto& n : navList) {
                if (n.sat == obs.sat && n.sys == obs.sys) {
                    double dt = timediff(obs.time, n.toe);
                    if (fabs(dt) < fabs(min_dt)) { min_dt = dt; nav = &n; }
                }
            }
            if (!nav || fabs(min_dt) > 7200.0) continue;

            // 2. 发射时刻 (扣除伪距传播时间)
            GPSTime t_tx = obs.time;
            t_tx.sec -= obs.P1 / CLIGHT;

            // 3. 卫星位置与钟差
            double rs[3], dts;
            satPos(t_tx, *nav, rs, &dts);

            // 4. 地球自转改正
            double r_dist = sqrt(pow(rs[0]-x[0],2) + pow(rs[1]-x[1],2) + pow(rs[2]-x[2],2));
            double tau = r_dist / CLIGHT;
            double theta = OMEGA_E * tau;
            double rx = rs[0]*cos(theta) + rs[1]*sin(theta);
            double ry = -rs[0]*sin(theta) + rs[1]*cos(theta);
            rs[0] = rx; rs[1] = ry; // 更新卫星坐标

            // 5. 几何距离与视线向量
            double r = sqrt(pow(rs[0]-x[0],2) + pow(rs[1]-x[1],2) + pow(rs[2]-x[2],2));
            double ex = (x[0] - rs[0]) / r;
            double ey = (x[1] - rs[1]) / r;
            double ez = (x[2] - rs[2]) / r;

            // 6. 高度角 (用于定权)
            double azel[2] = {0, PI/2};
            if (sqrt(x[0]*x[0]+x[1]*x[1]) > 1000.0) sat_azel(x, rs, azel);
            double el = azel[1];
            if (el < 10.0 * PI / 180.0) continue; // 截止高度角 10度

            // 7. 对流层改正 (Hopfield / Saastamoinen)
            double trop = 0.0, trop_var = 0.0;
            if (iter > 0) trop_model_prec(obs.time, x, azel, &trop, &trop_var);

            // 8. 组建方程
            // 观测值: 无电离层组合 IF
            double P_IF = ALPHA * obs.P1 - BETA * obs.P2;
            // 计算值: 几何距离 + 接收机钟差 - 卫星钟差 + 对流层
            double P_comp = r + *dtr - CLIGHT * dts + trop;

            double res = P_IF - P_comp;

            H.push_back(ex); H.push_back(ey); H.push_back(ez); H.push_back(1.0);
            v.push_back(res);

            // 简单定权模型: sin(el)^2
            W.push_back(sin(el)*sin(el));
            nSat++;
        }

        if (nSat < 4) return 0; // 卫星数不足

        // 9. 最小二乘解算
        int m = 4;
        vector<double> Q(m*m, 0.0);
        vector<double> b(m, 0.0);

        for(int i=0; i<nSat; i++) {
            double w = W[i];
            for(int j=0; j<4; j++) {
                for(int k=0; k<4; k++) Q[j*4+k] += H[i*4+j] * w * H[i*4+k];
                b[j] += H[i*4+j] * w * v[i];
            }
        }

        if (matinv(Q.data(), m) == -1) return 0; // 矩阵求逆失败

        double dx[4] = {0};
        matmul("NN", 4, 4, 1, 1.0, Q.data(), b.data(), 0.0, dx);

        // 更新状态
        x[0] += dx[0]; x[1] += dx[1]; x[2] += dx[2]; *dtr += dx[3];

        // 计算 PDOP
        *pdop = sqrt(Q[0] + Q[5] + Q[10]);

        // 收敛判据
        if (sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]) < 1e-4) return nSat;
    }
    return 0; // 未收敛
}

// ------------------------------------------------------------------------------------------------
// 核心模块 2: 多普勒测速 (Velocity Estimation)
// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
// 核心模块 2: 多普勒测速 (Velocity Estimation) - GPS双频IF组合版
// ------------------------------------------------------------------------------------------------
static int estVelocity(const vector<ObsData>& obsList, const vector<NavData>& navList,
                       const double* rr, double* vr, double* dtr_drift) {
    // vr: [vx, vy, vz], dtr_drift: [clk_drift]
    // 初始化
    double x[4] = {0, 0, 0, 0}; // 增量初值
    int nSat = 0;

    // 速度解算迭代
    for (int iter = 0; iter < 3; iter++) {
        vector<double> H, v, W;
        nSat = 0;

        for (const auto& obs : obsList) {
            // 1. 仅处理 GPS 系统
            if (obs.sys != 'G') continue;

            // 2. [关键修改] 检查是否具备双频多普勒数据
            // 如果要做严格的双频测速，必须 D1 和 D2 都有值
            if (obs.D1 == 0.0 || obs.D2 == 0.0) continue;

            // 3. 匹配星历
            const NavData* nav = nullptr;
            double min_dt = 1e9;
            for (const auto& n : navList) {
                if (n.sat == obs.sat && n.sys == obs.sys) {
                    double dt = timediff(obs.time, n.toe);
                    if (fabs(dt) < fabs(min_dt)) { min_dt = dt; nav = &n; }
                }
            }
            if (!nav) continue;

            // 4. 计算卫星状态
            GPSTime t_tx = obs.time;
            // 使用 IF 组合伪距估算传播时间，或默认值
            double dist_est = (obs.P1 > 0 && obs.P2 > 0) ? (ALPHA * obs.P1 - BETA * obs.P2) : obs.P1;
            if (dist_est > 0) t_tx.sec -= dist_est / CLIGHT;
            else t_tx.sec -= 0.070;

            double rs[3], vs[3], dts, dts_drift;
            satPosVel(t_tx, *nav, rs, vs, &dts, &dts_drift);

            // 5. 几何关系
            double r = sqrt(pow(rs[0]-rr[0],2) + pow(rs[1]-rr[1],2) + pow(rs[2]-rr[2],2));
            double e[3] = { (rs[0]-rr[0])/r, (rs[1]-rr[1])/r, (rs[2]-rr[2])/r };

            // 6. [关键修改] 计算双频无电离层组合伪距率 (Observed Range Rate)
            // 分别计算 L1 和 L2 的伪距率 (m/s)
            double rate_L1 = -obs.D1 * LAMBDA_L1;
            double rate_L2 = -obs.D2 * LAMBDA_L2;

            // IF 组合公式: rate_IF = alpha * rate_L1 - beta * rate_L2
            double obs_rate_IF = ALPHA * rate_L1 - BETA * rate_L2;

            // 7. 计算模型伪距率 (Computed Range Rate)
            double rate_sat_proj = e[0]*vs[0] + e[1]*vs[1] + e[2]*vs[2];

            // 预测值 = 卫地相对速度 + 钟漂差
            double rate_pred = rate_sat_proj - (e[0]*vr[0] + e[1]*vr[1] + e[2]*vr[2])
                               - CLIGHT * dts_drift + *dtr_drift;

            // 8. 残差与设计矩阵
            v.push_back(obs_rate_IF - rate_pred);

            // H 矩阵: [-ex, -ey, -ez, 1]
            H.push_back(-e[0]); H.push_back(-e[1]); H.push_back(-e[2]); H.push_back(1.0);

            // 定权 (简单等权或高度角权)
            W.push_back(1.0);
            nSat++;
        }

        if (nSat < 4) return 0;

        // 9. 最小二乘解算
        int m = 4;
        vector<double> Q(m*m, 0.0);
        vector<double> b(m, 0.0);

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
    cout << "MYGNSS10 SPP Processing (Position + Velocity)" << endl;
    cout << "Time Filter: 02:00:00 - 06:00:00 (GPST)" << endl;
    cout << "Ref Pos: " << fixed << setprecision(3) << refPos[0] << " " << refPos[1] << " " << refPos[2] << endl;
    cout << "--------------------------------------------------------------------------------------------------------" << endl;
    cout << "   GPST Time      Ns   PDOP    dE(m)    dN(m)    dU(m)   |   VE(m/s)  VN(m/s)  VU(m/s)" << endl;
    cout << "--------------------------------------------------------------------------------------------------------" << endl;
    cout << fixed << setprecision(3);

    for (const auto& epochObs : obsList) {
        if (epochObs.empty()) continue;
        GPSTime t_obs = epochObs[0].time;

        // --- 时间过滤逻辑 [新增] ---
        // 将 GPSTime 转换为历元数组 {year, month, day, hour, min, sec}
        double ep_check[6];
        time2epoch(t_obs, ep_check);

        // 计算当前小时数 (例如 2.5 表示 2:30)
        double current_hour = ep_check[3] + ep_check[4]/60.0 + ep_check[5]/3600.0;

        // 判断是否在 2:00 到 6:00 之间 (允许 1ms 的浮点数误差)
        // 范围: [2.0, 6.0]
        if (current_hour < 2.0 - 1e-9 || current_hour > 6.0 + 1e-9) {
            continue; // 不在范围内，跳过该历元
        }
        // ------------------------

        // 1. 定义状态变量
        double rr[3] = {0};     // 接收机位置
        double dtr = 0.0;       // 接收机钟差
        double pdop = 0.0;

        // 使用参考坐标作为先验初值
        if (refPos[0] != 0) { rr[0]=refPos[0]; rr[1]=refPos[1]; rr[2]=refPos[2]; }

        // 2. 执行定位
        int nSatPos = estPosition(epochObs, navList, rr, &dtr, &pdop);

        if (nSatPos >= 4) {
            // 定位成功，继续执行测速
            double vr[3] = {0};      // 接收机速度 [vx, vy, vz]
            double dtr_drift = 0.0;  // 接收机钟漂 (m/s)

            int nSatVel = estVelocity(epochObs, navList, rr, vr, &dtr_drift);

            // 3. 结果输出
            double ep[6];
            time2epoch(t_obs, ep);
            char time_str[32];
            sprintf(time_str, "%02.0f:%02.0f:%04.1f", ep[3], ep[4], ep[5]);

            // 计算 ENU 误差 (位置)
            double enu_pos[3] = {0};
            if (refPos[0] != 0) xyz2enu(rr, refPos, enu_pos);

            // 计算 ENU 速度
            double enu_vel[3] = {0};
            vel2enu(vr, rr, enu_vel);

            cout << time_str << "   "
                 << setw(2) << nSatPos << "   "
                 << setw(5) << pdop << "   "
                 << setw(8) << enu_pos[0] << " "
                 << setw(8) << enu_pos[1] << " "
                 << setw(8) << enu_pos[2] << "   |   "
                 << setw(7) << enu_vel[0] << "  "
                 << setw(7) << enu_vel[1] << "  "
                 << setw(7) << enu_vel[2] << endl;

        } else {
            // 定位失败
            double ep[6]; time2epoch(t_obs, ep);
            printf("%02.0f:%02.0f:%04.1f   SPP Failed (Ns=%d)\n", ep[3], ep[4], ep[5], nSatPos);
        }
    }
}
