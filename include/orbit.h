/*
* include/orbit.h - 卫星轨道与钟差计算声明
* 更新：增加 satPosVel 接口用于速度估计
*/

#ifndef MYGNSS10_ORBIT_H
#define MYGNSS10_ORBIT_H

#include "gtime.h" // 需要用到 GPSTime
#include "rinex.h" // 需要用到 NavData

// 定义光速 (供外部使用)
const double CLIGHT = 299792458.0;

// 地球自转角速度 (WGS84)
const double OMEGA_EARTH = 7.2921151467E-5;

/**
* 计算卫星位置和钟差 (仅位置，用于定位)
* @param tTx  信号发射时刻 (GPSTime)
* @param nav  对应的广播星历数据
* @param rs   [输出] 卫星位置 (ECEF坐标: x, y, z) 单位:米
* @param dts  [输出] 卫星钟差 (包含相对论效应) 单位:秒
*/
void satPos(GPSTime tTx, const NavData& nav, double* rs, double* dts);

/**
* 计算卫星位置、速度、钟差和钟漂 (用于测速)
* @param tTx       信号发射时刻 (GPSTime)
* @param nav       对应的广播星历数据
* @param rs        [输出] 卫星位置 (ECEF, 3x1) 单位:米
* @param vs        [输出] 卫星速度 (ECEF, 3x1) 单位:米/秒
* @param dts       [输出] 卫星钟差 单位:秒
* @param dts_drift [输出] 卫星钟漂 单位:秒/秒
*/
void satPosVel(GPSTime tTx, const NavData& nav, double* rs, double* vs, double* dts, double* dts_drift);

#endif // MYGNSS10_ORBIT_H