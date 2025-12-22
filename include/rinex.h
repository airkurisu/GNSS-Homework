/*
* rinex.h - RINEX 3.04 数据结构与读取接口声明
 * 更新：增加多普勒观测值 (D1, D2) 支持，用于速度解算
 */

#ifndef MYGNSS10_RINEX_H
#define MYGNSS10_RINEX_H

#include <vector>
#include <string>
#include <iostream>
#include "gtime.h"

// --- 导航星历结构体 (NavData) ---
// 对应广播星历文件 (.24p) 中的一条记录
struct NavData {
    int sat;        // 卫星号 (PRN)
    char sys;       // 卫星系统 (G:GPS, C:BDS, E:GAL, R:GLO)

    GPSTime toe;    // 星历参考时刻 (Time of Ephemeris)
    GPSTime toc;    // 钟差参考时刻 (Time of Clock)

    // 卫星钟差参数 (对应文件第1行)
    double a0, a1, a2;

    // --- 广播轨道参数 ---
    double crs, delta_n, M0;
    double cuc, e, cus, sqrtA;
    double cic, Omega0, cis;
    double i0, crc, omega, OmegaDot;
    double IDot;

    // --- 重要修正参数 ---
    double tgd;     // 群延迟 (Total Group Delay)
};

// --- 观测数据结构体 (ObsData) ---
// 对应观测文件 (.24o) 中某一历元下的一颗卫星
struct ObsData {
    GPSTime time;   // 观测时刻 (接收机时间)
    int sat;        // 卫星号
    char sys;       // 卫星系统 (G/C/E/R)

    // --- 伪距观测值 (单位: 米) ---
    // P1: L1/B1/E1/G1 频段
    // P2: L2/B3/E5a/G2 频段
    double P1;
    double P2;

    // --- 多普勒观测值 (单位: Hz) ---
    // 用于接收机速度估计
    double D1;
    double D2;
};

// --- 函数声明 ---

std::vector<NavData> readNavFile(const std::string& filename);
std::vector<std::vector<ObsData>> readObsFile(const std::string& filename, double* approxPos);

#endif // GNSS_SPP_RINEX_H