/*
 * src/rinex.cpp - RINEX Reader with Doppler Support
 * 更新说明:
 * 1. 保持了对 GPS/BDS/GAL/GLO 时间系统的归化处理。
 * 2. 增加了读取 'D' (Doppler) 类型观测值的逻辑。
 */

#include "rinex.h"
#include "gtime.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <vector>
#include <string>

using namespace std;

// 辅助函数：字符串转数字
static double str2num(const string& line, int start, int len) {
    if (start >= line.length()) return 0.0;
    string sub = line.substr(start, min(len, (int)line.length() - start));
    // 快速检查是否为空白
    if (sub.find_first_not_of(" \t\r\n") == string::npos) return 0.0;

    // 处理科学计数法中的 'D' 或 'd' (RINEX 格式遗留)
    size_t dPos = sub.find('D'); if (dPos != string::npos) sub[dPos] = 'E';
    dPos = sub.find('d'); if (dPos != string::npos) sub[dPos] = 'E';

    try { return stod(sub); } catch (...) { return 0.0; }
}

// 1. 读取导航电文 (Nav File)
vector<NavData> readNavFile(const string& filename) {
    vector<NavData> navList;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "[Error] Cannot open nav file: " << filename << endl;
        return navList;
    }

    string line;
    // 跳过文件头
    while (getline(file, line)) if (line.find("END OF HEADER") != string::npos) break;

    while (getline(file, line)) {
        if (line.length() < 20) continue;

        NavData nav = {0};
        string sysStr = line.substr(0, 1);
        nav.sys = sysStr[0];
        try { nav.sat = stoi(line.substr(1, 2)); } catch (...) { continue; }

        // 统一卫星编号
        if (nav.sys == 'C') nav.sat += 200;
        else if (nav.sys == 'E') nav.sat += 300;

        // 解析时间 (PRN epoch svclk)
        double ep[6];
        ep[0] = str2num(line, 4, 4); ep[1] = str2num(line, 9, 2); ep[2] = str2num(line, 12, 2);
        ep[3] = str2num(line, 15, 2); ep[4] = str2num(line, 18, 2); ep[5] = str2num(line, 21, 2);

        nav.toc = epoch2time(ep);

        // --- 时间系统归化 (统一到 GPST) ---
        if (nav.sys == 'R') {
            // GLONASS UTC -> GPST (跳秒+18s, 近似处理)
            nav.toc = utc2gpst(nav.toc);
        }
        else if (nav.sys == 'C') {
            // BDS BDT -> GPST (+14s)
            nav.toc.sec += 14.0;
            if (nav.toc.sec >= 604800.0) { nav.toc.sec -= 604800.0; nav.toc.week++; }
        }

        // 读取不同系统的轨道参数
        if (nav.sys == 'R') {
            // GLONASS: TauN, GammaN
            nav.a0  = str2num(line, 23, 19);
            nav.a1  = str2num(line, 42, 19);
            nav.toe = nav.toc;

            if (!getline(file, line)) break;
            nav.M0     = str2num(line, 4, 19) * 1000.0;  // X
            nav.Omega0 = str2num(line, 23, 19) * 1000.0; // Vx
            nav.cuc    = str2num(line, 42, 19) * 1000.0; // Ax

            if (!getline(file, line)) break;
            nav.e      = str2num(line, 4, 19) * 1000.0;  // Y
            nav.i0     = str2num(line, 23, 19) * 1000.0; // Vy
            nav.cus    = str2num(line, 42, 19) * 1000.0; // Ay
            nav.a2     = str2num(line, 61, 19);          // FCN

            if (!getline(file, line)) break;
            nav.sqrtA  = str2num(line, 4, 19) * 1000.0;  // Z
            nav.omega  = str2num(line, 23, 19) * 1000.0; // Vz
            nav.crc    = str2num(line, 42, 19) * 1000.0; // Az

        } else {
            // GPS/BDS/GAL: Kepler
            nav.a0 = str2num(line, 23, 19); nav.a1 = str2num(line, 42, 19); nav.a2 = str2num(line, 61, 19);

            if (!getline(file, line)) break;
            nav.crs = str2num(line, 23, 19); nav.delta_n = str2num(line, 42, 19); nav.M0 = str2num(line, 61, 19);

            if (!getline(file, line)) break;
            nav.cuc = str2num(line, 4, 19); nav.e = str2num(line, 23, 19); nav.cus = str2num(line, 42, 19); nav.sqrtA = str2num(line, 61, 19);

            if (!getline(file, line)) break;
            nav.toe.sec = str2num(line, 4, 19);
            nav.toe.week = nav.toc.week;
            // 简单处理周跨越
            double dt = nav.toe.sec - nav.toc.sec;
            if (dt > 302400.0) nav.toe.week -= 1; else if (dt < -302400.0) nav.toe.week += 1;

            nav.cic = str2num(line, 23, 19); nav.Omega0 = str2num(line, 42, 19); nav.cis = str2num(line, 61, 19);

            if (!getline(file, line)) break;
            nav.i0 = str2num(line, 4, 19); nav.crc = str2num(line, 23, 19); nav.omega = str2num(line, 42, 19); nav.OmegaDot = str2num(line, 61, 19);

            if (!getline(file, line)) break;
            nav.IDot = str2num(line, 4, 19);

            if (!getline(file, line)) break;
            // TGD 位置可能不同
            if (nav.sys == 'G') nav.tgd = str2num(line, 42, 19);
            else if (nav.sys == 'C') nav.tgd = str2num(line, 4, 19); // TGD1 (B1/B3)
            else if (nav.sys == 'E') nav.tgd = str2num(line, 42, 19); // BGD E1/E5a
        }
        navList.push_back(nav);
    }
    file.close();
    cout << "Read " << navList.size() << " nav records." << endl;
    return navList;
}

// 2. 读取观测文件 (Obs File)
vector<vector<ObsData>> readObsFile(const string& filename, double* approxPos) {
    vector<vector<ObsData>> allEpochs;
    ifstream file(filename);
    if (!file.is_open()) return allEpochs;

    string line;
    // 映射表: [系统][观测类型字符串] -> 数据列索引
    map<char, map<string, int>> sysObsTypeMap;
    // 记录每个系统总共有多少种观测类型
    map<char, int> sysObsCount;

    // --- 解析 Header ---
    while (getline(file, line)) {
        // 读取近似坐标 (Approx Position)
        if (line.find("APPROX POSITION XYZ") != string::npos) {
            approxPos[0] = str2num(line, 0, 14);
            approxPos[1] = str2num(line, 14, 14);
            approxPos[2] = str2num(line, 28, 14);
        }
        // 读取观测类型 (OBS TYPES)
        if (line.find("SYS / # / OBS TYPES") != string::npos) {
            char sys = line[0]; // G/C/E/R
            if (sys == 'G' || sys == 'C' || sys == 'E' || sys == 'R') {
                int nObs = stoi(line.substr(4, 2));
                string typesStr = line.substr(7, 53); // RINEX 3.04 格式

                // 处理续行 (若 nObs 很大，可能会跨行，此处暂简化处理单行或标准情况)
                // 在严谨的 RINEX 解析中需要循环读取直到 nObs 读完，这里假设一行够用或项目数据简单

                int currentIdx = sysObsCount[sys];
                for (int i = 0; i < typesStr.length(); i += 4) {
                    if (currentIdx >= nObs) break;
                    string type = typesStr.substr(i, 3);
                    if (isalpha(type[0])) {
                        sysObsTypeMap[sys][type] = currentIdx;
                        currentIdx++;
                    }
                }
                sysObsCount[sys] = currentIdx;
            }
        }
        if (line.find("END OF HEADER") != string::npos) break;
    }

    // --- 解析 Body (Epochs) ---
    while (getline(file, line)) {
        // 查找历元开始标记 '>'
        if (line.empty() || line[0] != '>') continue;

        // 解析历元时间
        double ep[6];
        ep[0] = str2num(line, 2, 4); ep[1] = str2num(line, 7, 2); ep[2] = str2num(line, 10, 2);
        ep[3] = str2num(line, 13, 2); ep[4] = str2num(line, 16, 2); ep[5] = str2num(line, 19, 11);
        GPSTime t = epoch2time(ep);

        int nSat = 0;
        try { nSat = stoi(line.substr(32, 3)); } catch (...) { nSat = 0; }

        vector<ObsData> epochObs;

        for (int i = 0; i < nSat; i++) {
            if (!getline(file, line)) break;

            char sys = line[0];
            if (sysObsTypeMap.find(sys) == sysObsTypeMap.end()) continue; // 跳过未知系统

            ObsData obs = {0};
            obs.time = t;
            obs.sys = sys;

            int prn = 0;
            try { prn = stoi(line.substr(1, 2)); } catch (...) { continue; }

            // 统一卫星号
            if (sys == 'C') obs.sat = prn + 200;
            else if (sys == 'E') obs.sat = prn + 300;
            else if (sys == 'R') obs.sat = prn; // GLONASS 保持原号或按需+偏移
            else obs.sat = prn;

            auto& typeMap = sysObsTypeMap[sys];

            // --- 提取观测值 (P:伪距, D:多普勒) ---
            // 策略: 优先读取特定码, 若无则尝试备选

            // GPS (L1, L2)
            if (sys == 'G') {
                // P1
                if (typeMap.count("C1C")) obs.P1 = str2num(line, 3 + typeMap["C1C"] * 16, 14);
                else if (typeMap.count("C1W")) obs.P1 = str2num(line, 3 + typeMap["C1W"] * 16, 14);
                // P2
                if (typeMap.count("C2W")) obs.P2 = str2num(line, 3 + typeMap["C2W"] * 16, 14);
                else if (typeMap.count("C2L")) obs.P2 = str2num(line, 3 + typeMap["C2L"] * 16, 14);

                // D1
                if (typeMap.count("D1C")) obs.D1 = str2num(line, 3 + typeMap["D1C"] * 16, 14);
                else if (typeMap.count("D1W")) obs.D1 = str2num(line, 3 + typeMap["D1W"] * 16, 14);
                // D2
                if (typeMap.count("D2W")) obs.D2 = str2num(line, 3 + typeMap["D2W"] * 16, 14);
                else if (typeMap.count("D2L")) obs.D2 = str2num(line, 3 + typeMap["D2L"] * 16, 14);
            }
            // BDS (B1, B3) -> 对应 RINEX 3 中的 C2I, C6I
            else if (sys == 'C') {
                // P1 (B1)
                if (typeMap.count("C2I")) obs.P1 = str2num(line, 3 + typeMap["C2I"] * 16, 14);
                // P2 (B3)
                if (typeMap.count("C6I")) obs.P2 = str2num(line, 3 + typeMap["C6I"] * 16, 14);

                // D1 (B1)
                if (typeMap.count("D2I")) obs.D1 = str2num(line, 3 + typeMap["D2I"] * 16, 14);
                // D2 (B3)
                if (typeMap.count("D6I")) obs.D2 = str2num(line, 3 + typeMap["D6I"] * 16, 14);
            }
            // Galileo (E1, E5a)
            else if (sys == 'E') {
                if (typeMap.count("C1X")) obs.P1 = str2num(line, 3 + typeMap["C1X"] * 16, 14);
                if (typeMap.count("C5X")) obs.P2 = str2num(line, 3 + typeMap["C5X"] * 16, 14);

                if (typeMap.count("D1X")) obs.D1 = str2num(line, 3 + typeMap["D1X"] * 16, 14);
                if (typeMap.count("D5X")) obs.D2 = str2num(line, 3 + typeMap["D5X"] * 16, 14);
            }
            // GLONASS (G1, G2)
            else if (sys == 'R') {
                if (typeMap.count("C1C")) obs.P1 = str2num(line, 3 + typeMap["C1C"] * 16, 14);
                else if (typeMap.count("C1P")) obs.P1 = str2num(line, 3 + typeMap["C1P"] * 16, 14);

                if (typeMap.count("C2C")) obs.P2 = str2num(line, 3 + typeMap["C2C"] * 16, 14);
                else if (typeMap.count("C2P")) obs.P2 = str2num(line, 3 + typeMap["C2P"] * 16, 14);

                if (typeMap.count("D1C")) obs.D1 = str2num(line, 3 + typeMap["D1C"] * 16, 14);
                else if (typeMap.count("D1P")) obs.D1 = str2num(line, 3 + typeMap["D1P"] * 16, 14);

                if (typeMap.count("D2C")) obs.D2 = str2num(line, 3 + typeMap["D2C"] * 16, 14);
                else if (typeMap.count("D2P")) obs.D2 = str2num(line, 3 + typeMap["D2P"] * 16, 14);
            }

            // 只有当至少有一个伪距或一个多普勒值有效时才存储
            // (通常解算至少需要伪距，如果只做测速可能只需要D，但这里做SPP+Vel，假设需要伪距)
            if (obs.P1 > 0.1 || obs.P2 > 0.1) {
                epochObs.push_back(obs);
            }
        }
        if (!epochObs.empty()) allEpochs.push_back(epochObs);
    }
    file.close();
    cout << "Read " << allEpochs.size() << " epochs." << endl;
    return allEpochs;
}