/*
* include/rinex.h - RINEX 3.04 数据结构
* 修复记录:
* 1. 保持结构体对齐，确保 NavData 包含 freq_num
*/

#ifndef MYGNSS10_RINEX_H
#define MYGNSS10_RINEX_H

#include <vector>
#include <string>
#include <map>
#include "gtime.h"

// --- 导航星历结构体 ---
struct NavData {
    int sat;        // PRN
    char sys;       // System (G/C/E/R)
    GPSTime toe;    // Time of Ephemeris
    GPSTime toc;    // Time of Clock

    // 钟差参数
    double a0, a1, a2;

    // 轨道参数 (兼容 Kepler 和 GLONASS State Vector)
    // 变量名复用以保持内存紧凑
    double crs, delta_n, M0;   // GLO: X, Vx, Ax
    double cuc, e, cus, sqrtA; // GLO: Y, Vy, Ay, FreqNum (via helper)
    double cic, Omega0, cis;   // GLO: Z, Vz, Az
    double i0, crc, omega, OmegaDot;
    double IDot;

    double tgd;     // TGD (BDS/GPS)
    int freq_num;   // GLONASS 频率号 (关键: 用于计算波长)
};

// --- 观测数据结构体 ---
struct ObsData {
    GPSTime time;
    int sat;
    char sys;

    // 核心观测值
    double P1, P2;  // Pseudorange (m)
    double D1, D2;  // Doppler (Hz)

    // 观测码类型 (调试用)
    std::string codeP1, codeP2;
};

// --- 函数声明 ---
std::vector<NavData> readNavFile(const std::string& filename);
std::vector<std::vector<ObsData>> readObsFile(const std::string& filename, double* approxPos);

#endif // MYGNSS10_RINEX_H