/*
* include/tropo.h - 对流层延迟改正模块声明
 * 模型: Improved Hopfield (ZTD) + NMF (Mapping Function)
 */

#ifndef MYGNSS10_TROPO_H
#define MYGNSS10_TROPO_H

#include "gtime.h" // 需要用到时间结构体

// --- 常量定义 ---
// 标准大气参数 (Standard Atmosphere)
const double PR_STD = 1013.25;      // 标准气压 (hPa)
const double TR_STD = 288.15;       // 标准温度 (K)
const double HR_STD = 50.0;         // 标准相对湿度 (%)

// --- 函数声明 ---

/**
 * 对流层模型主函数 (Improved Hopfield + NMF)
 * @param time   当前观测时间 (用于NMF计算年积日)
 * @param pos    测站坐标 {纬度(rad), 经度(rad), 高程(m)}
 * @param azel   卫星方位角、高度角 {az, el} (rad)
 * @param trop   [输出] 对流层总延迟 (m)
 * @param var    [输出] 对流层误差方差 (m^2)
 * @return       1: 成功, 0: 失败(高度角过低)
 */
int trop_model_prec(GPSTime time, const double *pos, const double *azel,
                    double *trop, double *var);

#endif // GNSS_SPP_TROPO_H