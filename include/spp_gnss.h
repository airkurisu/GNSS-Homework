/*
* include/spp_gnss.h - 多系统GNSS SPP定位与测速声明
*/

#ifndef MYGNSS10_SPP_GNSS_H
#define MYGNSS10_SPP_GNSS_H

#include <vector>
#include "rinex.h"

/**
 * 执行 多系统 GNSS SPP 定位与测速
 * 支持 GPS, BDS, Galileo 等 (需根据观测文件实际包含的系统)
 * 自动处理系统间偏差 (ISB)
 *
 * @param obsList  观测数据
 * @param navList  星历数据
 * @param refPos   参考坐标
 */
void spp_gnss_process(const std::vector<std::vector<ObsData>>& obsList,
                      const std::vector<NavData>& navList,
                      const double* refPos);

#endif // MYGNSS10_SPP_GNSS_H