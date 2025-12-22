/*
* include/spp.h - SPP定位与测速核心声明
* 项目: MYGNSS10
*/

#ifndef MYGNSS10_SPP_H
#define MYGNSS10_SPP_H

#include <vector>
#include "rinex.h" // 需要 ObsData, NavData

// 定义定位结果结构体（可选，用于存储单历元结果）
struct SPPResult {
    GPSTime time;
    double pos[3];  // ECEF x, y, z (m)
    double vel[3];  // ECEF vx, vy, vz (m/s)
    double dtr;     // 接收机钟差 (s)
    double dtr_rate;// 接收机钟漂 (s/s)
    int nSat;
    double pdop;
};

/**
 * 执行 SPP 定位与测速主流程
 * @param obsList  所有历元的观测数据
 * @param navList  导航星历数据
 * @param refPos   参考坐标 (ECEF XYZ, 单位:米)，用于计算 ENU 误差，若无则传 {0,0,0}
 */
void spp_process(const std::vector<std::vector<ObsData>>& obsList,
                 const std::vector<NavData>& navList,
                 const double* refPos);

#endif // MYGNSS10_SPP_H