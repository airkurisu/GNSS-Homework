/*
 * src/tropo.cpp - 对流层延迟改正具体实现
 * * 算法参考文献:
 * 1. Hopfield Improved Model: Goad and Goodman (1974)
 * 2. NMF (Niell Mapping Function): Niell (1996)
 */

#include "tropo.h"
#include <cmath>
#include <iostream>

using namespace std;

const double PI = 3.1415926535897932;

// --- 内部辅助函数: 获取年积日 (Day of Year, DOY) ---
// 依赖 gtime.h 中的 time2epoch
static double time2doy(GPSTime t) {
    double ep[6];
    time2epoch(t, ep);
    // 简单计算 DOY (不考虑复杂的闰年平年细微差异，对NMF精度足够)
    // NMF 公式的 DOY 是从 0 还是 1 开始差异极小，这里采用标准算法
    // 简易算法：(time - year_start)
    double epoch_first[6] = {ep[0], 1, 1, 0, 0, 0};
    GPSTime t_first = epoch2time(epoch_first);
    double dt = timediff(t, t_first);
    return dt / 86400.0 + 1.0;
}

// --- 1. 标准大气模型 (Standard Atmosphere) ---
// 当没有实测气象参数时，根据高程估算 P, T, H
static void get_met_standard(double h, double *pres, double *temp, double *e) {
    // h: 测站高程 (m)
    // pres: 气压 (hPa), temp: 温度 (K), e: 水汽压 (hPa)

    // 限制高度范围 (-100m ~ 10000m)
    if (h < -100.0) h = -100.0;
    if (h > 10000.0) h = 10000.0;

    // 标准大气递减率
    *pres = PR_STD * pow(1.0 - 0.0000226 * h, 5.225);
    *temp = TR_STD - 0.0065 * h;
    double humidity = HR_STD / 100.0; // 相对湿度 0.5

    // 计算饱和水汽压 (Magnus 公式)
    // es = 6.11 * 10 ^ (7.5 * Tc / (237.3 + Tc)), Tc = T - 273.15
    double Tc = *temp - 273.15;
    double es = 6.11 * pow(10.0, (7.5 * Tc) / (237.3 + Tc));

    // 实际水汽压
    *e = humidity * es;
}

// --- 2. Hopfield 改进模型 (计算天顶延迟 ZHD, ZWD) ---
static void trop_hopfield_zenith(double P, double T, double e,
                                 double *zhd, double *zwd)
{
    // P(hPa), T(K), e(hPa)

    // 干分量折射率参数
    // Nd = 77.64 * (P/T)
    // 湿分量折射率参数
    // Nw = -12.96 * (e/T) + 3.718e5 * (e/T^2)

    // 改进 Hopfield 的对流层顶高度 (m)
    double h_dry = 40136.0 + 148.72 * (T - 273.16);
    double h_wet = 11000.0; // 湿分量高度通常固定为 11km 或 12km

    // 计算天顶延迟 (Zenith Delay)
    // 公式: D_z = (10^-6 / 5) * N_0 * h_eff

    double N_d0 = 77.64 * (P / T);
    double N_w0 = -12.96 * (e / T) + 3.718e5 * (e / (T * T));

    *zhd = 1.552e-5 * (P / T) * h_dry; // 化简后的公式
    *zwd = 1.552e-5 * (N_w0 / 77.64) * h_wet; // 类似的比例关系

    // 或者使用更严谨的积分化简形式:
    // ZHD = 1e-6 * N_d0 * h_dry / 5.0;
    // ZWD = 1e-6 * N_w0 * h_wet / 5.0;
    // 上面的 1.552e-5 其实就是 77.6 * 1e-6 / 5 * h_factor...
    // 为了严谨，我们用标准展开式：

    *zhd = 1.0e-6 / 5.0 * N_d0 * h_dry;
    *zwd = 1.0e-6 / 5.0 * N_w0 * h_wet;
}

// --- 3. NMF 投影函数 (Niell Mapping Function) ---

// NMF 连分式插值函数
static double nmf_func(double el, double a, double b, double c) {
    double sin_el = sin(el);
    double top = 1.0 + a / (1.0 + b / (1.0 + c));
    double bot = sin_el + a / (sin_el + b / (sin_el + c));
    return top / bot;
}

// NMF 系数表 (按纬度 15, 30, 45, 60, 75)
static const double COEF_DRY[5][3][2] = {
    // lat 15: {a_avg, a_amp}, {b_avg, b_amp}, {c_avg, c_amp}
    {{1.2769934e-3, 0.0}, {2.9153695e-3, 0.0}, {62.610505e-3, 0.0}},
    // lat 30:
    {{1.2683230e-3, 1.2709626e-5}, {2.9152299e-3, 2.1414979e-5}, {62.837393e-3, 9.0128400e-5}},
    // lat 45:
    {{1.2465397e-3, 2.6523662e-5}, {2.9288445e-3, 3.0160779e-5}, {63.721774e-3, 4.3497037e-5}},
    // lat 60:
    {{1.2196049e-3, 3.4000452e-5}, {2.9022565e-3, 7.2543778e-5}, {63.824265e-3, 8.4795348e-5}},
    // lat 75:
    {{1.2045996e-3, 4.1202191e-5}, {2.9024912e-3, 11.723375e-5}, {64.258455e-3, 17.037206e-5}}
};

static const double COEF_WET[5][3] = {
    // a, b, c (无季节变化)
    {5.8021897e-4, 1.4275268e-3, 4.3472964e-2}, // 15
    {5.6794847e-4, 1.5138625e-3, 4.6729510e-2}, // 30
    {5.8118019e-4, 1.4572752e-3, 4.3908931e-2}, // 45
    {5.9727542e-4, 1.5007428e-3, 4.4626982e-2}, // 60
    {6.1641693e-4, 1.7599082e-3, 5.4736038e-2}  // 75
};

static void trop_nmf(GPSTime time, const double *pos, double el,
                     double *map_dry, double *map_wet)
{
    // pos: {lat, lon, h} (rad, rad, m)
    double lat = fabs(pos[0] * 180.0 / PI); // 纬度转为度，取绝对值
    double h = pos[2];
    double doy = time2doy(time); // 年积日

    // 1. 确定纬度插值索引
    int i;
    if      (lat <= 15.0) i = 0;
    else if (lat <= 30.0) i = 1;
    else if (lat <= 45.0) i = 2;
    else if (lat <= 60.0) i = 3;
    else                  i = 4; // lat > 75 实际上通常用 75 的值或外推，这里简化为用 i=4

    // 2. 计算干分量系数 (Dry Coefficients)
    // 需要对纬度进行线性插值
    double a_d, b_d, c_d;
    double zeros[3] = {0}; // 临时用

    // 为了简化代码，这里实现简单的线性插值逻辑
    // 实际 NMF 逻辑：如果 lat 在 15-75 之间，在两个节点间插值
    // 公式: val = val_i + (val_{i+1} - val_i) * (lat - lat_i) / (lat_{i+1} - lat_i)

    double coef_d[3]; // 最终的 a, b, c

    // 定义节点纬度
    double lats[] = {15.0, 30.0, 45.0, 60.0, 75.0};

    for (int k = 0; k < 3; k++) { // 遍历 a, b, c
        double val_avg, val_amp;

        if (lat <= 15.0) {
            val_avg = COEF_DRY[0][k][0];
            val_amp = COEF_DRY[0][k][1];
        } else if (lat >= 75.0) {
            val_avg = COEF_DRY[4][k][0];
            val_amp = COEF_DRY[4][k][1];
        } else {
            // 线性插值
            int idx = (int)((lat - 15.0) / 15.0); // 0..3
            double ratio = (lat - lats[idx]) / 15.0;

            double avg1 = COEF_DRY[idx][k][0];
            double avg2 = COEF_DRY[idx+1][k][0];
            double amp1 = COEF_DRY[idx][k][1];
            double amp2 = COEF_DRY[idx+1][k][1];

            val_avg = avg1 + (avg2 - avg1) * ratio;
            val_amp = amp1 + (amp2 - amp1) * ratio;
        }


        double sign = (pos[0] >= 0) ? 1.0 : -1.0;
        // 原始 NMF 公式是 P = P_avg - P_amp * cos(...)
        // 这里的实现我们遵循 RTKLIB/Gipsy 的标准:
        // 北半球: - amp * cos; 南半球: - amp * cos( ... + PI) = + amp * cos

        coef_d[k] = val_avg - sign * val_amp * cos(2.0 * PI * (doy - 28.0) / 365.25);
    }

    // 3. 计算干分量投影 (Height Correction)
    // NMF 干分量有高度修正项: dm = (1/sin(e) - f(e, a_ht, b_ht, c_ht)) * H / 1000
    // 高度修正系数 (固定值)
    const double a_ht = 2.53e-5;
    const double b_ht = 5.49e-3;
    const double c_ht = 1.14e-3;

    double m_dry_0 = nmf_func(el, coef_d[0], coef_d[1], coef_d[2]);
    double m_ht    = nmf_func(el, a_ht, b_ht, c_ht);

    *map_dry = m_dry_0 + (1.0 / sin(el) - m_ht) * (h / 1000.0);

    // 4. 计算湿分量系数 (Wet Coefficients)
    // 同样需要对纬度插值
    double coef_w[3];
    for (int k = 0; k < 3; k++) {
        if (lat <= 15.0) {
            coef_w[k] = COEF_WET[0][k];
        } else if (lat >= 75.0) {
            coef_w[k] = COEF_WET[4][k];
        } else {
            int idx = (int)((lat - 15.0) / 15.0);
            double ratio = (lat - lats[idx]) / 15.0;
            double val1 = COEF_WET[idx][k];
            double val2 = COEF_WET[idx+1][k];
            coef_w[k] = val1 + (val2 - val1) * ratio;
        }
    }

    *map_wet = nmf_func(el, coef_w[0], coef_w[1], coef_w[2]);
}


// --- 对流层模型主入口 ---
int trop_model_prec(GPSTime time, const double *pos, const double *azel,
                    double *trop, double *var)
{
    // 1. 检查高度角
    double el = azel[1];
    if (el < 3.0 * PI / 180.0) { // 极低高度角不计算
        *trop = 0.0;
        if (var) *var = 0.0;
        return 0;
    }

    // 2. 获取气象参数 (使用标准大气模型估算)
    // 这里的 pos[2] 是椭球高，严格来说应该用正高/海拔高，但在 SPP 中混用误差可忽略
    double pres, temp, e;
    get_met_standard(pos[2], &pres, &temp, &e);

    // 3. 计算天顶延迟 (Improved Hopfield)
    double zhd, zwd;
    trop_hopfield_zenith(pres, temp, e, &zhd, &zwd);

    // 4. 计算投影函数 (NMF)
    double map_dry, map_wet;
    trop_nmf(time, pos, el, &map_dry, &map_wet);

    // 5. 组合总延迟
    // Total = ZHD * M_dry + ZWD * M_wet
    *trop = zhd * map_dry + zwd * map_wet;

    // 6. 估算方差 (可选)
    // 简单给一个定值，约 0.3m 的 ZTD 误差映射到视线方向
    if (var) *var = 0.3 * 0.3;

    return 1;
}