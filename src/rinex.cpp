/*
 * src/rinex.cpp - RINEX 读取 (完整修复版)
 * 修正记录:
 * 1. [CRITICAL] 修复 GLONASS 时间系统问题 (UTC -> GPST)。
 * 2. 修复 GLONASS 轨道单位 (km -> m)。
 * 3. 包含完整的频率号读取逻辑。
 */

#include "rinex.h"
#include "gtime.h" 
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <map>
#include <iostream>
#include <vector>

using namespace std;

// --- 内部辅助结构与变量 ---

// 用于记录信号列索引
struct ColIdx {
    int p1=-1, p2=-1, d1=-1, d2=-1;
    string c1, c2;
};

// 信号优先级列表
// 注意：对于 GLONASS，通常使用 C1C/C2C 或 C1P/C2P
static map<char, vector<string>> prio_P1 = {
    {'G', {"C1C", "C1W", "C1X"}}, {'E', {"C1X", "C1C"}},
    {'C', {"C2I", "C1P", "C1X"}}, {'R', {"C1C", "C1P"}} 
};
static map<char, vector<string>> prio_P2 = {
    {'G', {"C2W", "C2P", "C2X"}}, {'E', {"C5X", "C5Q"}},
    {'C', {"C6I", "C7I"}},        {'R', {"C2C", "C2P"}}
};

// 字符串转数字 (处理 1.23D-04 格式)
static double str2num(const string& s) {
    try {
        if (s.find_first_not_of(" \t\r\n") == string::npos) return 0.0;
        string temp = s;
        for (char &c : temp) {
            if (c == 'D' || c == 'd') c = 'E';
        }
        return stod(temp);
    } catch (...) {
        return 0.0;
    }
}

// --------------------------------------------------------------------------
// 1. 导航电文读取
// --------------------------------------------------------------------------
vector<NavData> readNavFile(const string& filename) {
    vector<NavData> navList;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open nav file " << filename << endl;
        return navList;
    }

    string line;
    // 跳过 Header
    while (getline(file, line)) {
        if (line.find("END OF HEADER") != string::npos) break;
    }

    NavData nav = {0};
    vector<string> buff;

    while (getline(file, line)) {
        if (line.empty()) continue;

        // 判断是否为新记录 (非空格开头且长度足够)
        // RINEX NAV 文件中，第一列是卫星号，后面紧接时间
        if (line[0] != ' ' && line.length() > 5) {
            // 保存上一条
            if (nav.sat != 0) navList.push_back(nav);

            nav = {0};
            nav.sys = line[0];
            
            // 容错处理：有些文件可能是 " 1" 或 "G01"
            string sat_str = line.substr(1, 2);
            try { 
                nav.sat = stoi(sat_str); 
            } catch(...) { 
                // 尝试处理 " 1" 这种格式
                if (sat_str[0] == ' ') sat_str[0] = '0';
                try { nav.sat = stoi(sat_str); } catch(...) { continue; }
            }

            // 解析时间
            try {
                int yr = stoi(line.substr(4, 4));
                int mon = stoi(line.substr(9, 2));
                int day = stoi(line.substr(12, 2));
                int hr = stoi(line.substr(15, 2));
                int min = stoi(line.substr(18, 2));
                double sec = str2num(line.substr(21, 4));

                double toc_ep[6] = {(double)yr, (double)mon, (double)day, (double)hr, (double)min, sec};
                nav.toc = epoch2time(toc_ep);

                // =========================================================
                // [FIX] 关键修正: GLONASS 时间转换
                // GLONASS 导航电文时间 (TOC) 是 UTC，必须转为 GPST
                // 否则会导致轨道计算产生 ~18秒 (即 ~60km) 的误差
                // =========================================================
                if (nav.sys == 'R') {
                    nav.toc = utc2gpst(nav.toc);
                }
                // =========================================================

            } catch(...) { continue; }

            // 钟差参数
            nav.a0 = str2num(line.substr(23, 19)); // GPS: Clock Bias, GLO: -TauN
            nav.a1 = str2num(line.substr(42, 19)); // GPS: Drift,      GLO: +GammaN
            nav.a2 = str2num(line.substr(61, 19)); // GPS: Drift Rate, GLO: Message Frame Time
            buff.clear();
        } else {
            // 数据行
            buff.push_back(line);

            if (nav.sys == 'R') {
                // GLONASS (3行数据)
                // 必须将单位从 km 转为 m
                if (buff.size() >= 3) {
                    // Line 1 (Buffer 0): X
                    nav.M0     = str2num(buff[0].substr( 4, 19)) * 1000.0; // Pos X
                    nav.Omega0 = str2num(buff[0].substr(23, 19)) * 1000.0; // Vel X
                    nav.cuc    = str2num(buff[0].substr(42, 19)) * 1000.0; // Acc X

                    // Line 2 (Buffer 1): Y
                    nav.e      = str2num(buff[1].substr( 4, 19)) * 1000.0; // Pos Y
                    nav.i0     = str2num(buff[1].substr(23, 19)) * 1000.0; // Vel Y
                    nav.cus    = str2num(buff[1].substr(42, 19)) * 1000.0; // Acc Y

                    // 读取 GLONASS 频率号 (RINEX 2 标准位置: Line 2, Field 4)
                    if (buff[1].length() >= 61) {
                        nav.freq_num = (int)str2num(buff[1].substr(61, 19));
                    }

                    // Line 3 (Buffer 2): Z
                    nav.sqrtA  = str2num(buff[2].substr( 4, 19)) * 1000.0; // Pos Z
                    nav.omega  = str2num(buff[2].substr(23, 19)) * 1000.0; // Vel Z
                    nav.crc    = str2num(buff[2].substr(42, 19)) * 1000.0; // Acc Z

                    nav.toe = nav.toc; // GLONASS TOE = TOC
                }
            } else {
                // GPS/BDS/GAL (7行数据)
                if (buff.size() >= 7) {
                    nav.crs = str2num(buff[0].substr(23, 19));
                    nav.delta_n = str2num(buff[0].substr(42, 19));
                    nav.M0 = str2num(buff[0].substr(61, 19));

                    nav.cuc = str2num(buff[1].substr(4, 19));
                    nav.e = str2num(buff[1].substr(23, 19));
                    nav.cus = str2num(buff[1].substr(42, 19));
                    nav.sqrtA = str2num(buff[1].substr(61, 19));

                    nav.toe.sec = str2num(buff[2].substr(4, 19));
                    nav.toe.week = nav.toc.week; // 初步对齐到 TOC 周
                    nav.cic = str2num(buff[2].substr(23, 19));
                    nav.Omega0 = str2num(buff[2].substr(42, 19));
                    nav.cis = str2num(buff[2].substr(61, 19));

                    nav.i0 = str2num(buff[3].substr(4, 19));
                    nav.crc = str2num(buff[3].substr(23, 19));
                    nav.omega = str2num(buff[3].substr(42, 19));
                    nav.OmegaDot = str2num(buff[3].substr(61, 19));

                    nav.IDot = str2num(buff[4].substr(4, 19));
                    
                    // TGD 读取 (BDS 需要)
                    if (buff.size() >= 6) nav.tgd = str2num(buff[5].substr(42, 19));
                    
                    // 修正 TOE 周跨越 (如果 toe 与 toc 跨周)
                    double dt = nav.toe.sec - nav.toc.sec;
                    if (dt > 302400.0) nav.toe.week--;
                    else if (dt < -302400.0) nav.toe.week++;
                }
            }
        }
    }
    // 添加最后一个卫星
    if (nav.sat != 0) navList.push_back(nav);
    file.close();
    cout << "Loaded " << navList.size() << " navigation records." << endl;
    return navList;
}

// --------------------------------------------------------------------------
// 2. 观测文件读取
// --------------------------------------------------------------------------

// 辅助函数: 映射列号
static ColIdx map_signal_col(char sys, const vector<string>& types) {
    ColIdx idx;
    // 查找 P1 (优先频点)
    for (const string& target : prio_P1[sys]) {
        auto it = find(types.begin(), types.end(), target);
        if (it != types.end()) {
            idx.p1 = distance(types.begin(), it);
            idx.c1 = target;
            break;
        }
    }
    // Fallback P1: 如果没找到定义的优先级，找任意 C1/C2
    if (idx.p1 == -1) {
        for (size_t i=0; i<types.size(); i++) {
            if(types[i][0]!='C') continue;
            // 简单策略: BDS用C2I/C1x, 其它用C1x
            if((sys=='C' && (types[i][1]=='2'||types[i][1]=='1')) || (sys!='C' && types[i][1]=='1')) {
                idx.p1 = i; idx.c1 = types[i]; break;
            }
        }
    }
    // 查找 D1 (对应 P1)
    if (idx.p1 != -1) {
        string tD = idx.c1; tD[0] = 'D'; // 将 Cxx 改为 Dxx
        auto it = find(types.begin(), types.end(), tD);
        if (it != types.end()) idx.d1 = distance(types.begin(), it);
    }

    // 查找 P2
    for (const string& target : prio_P2[sys]) {
        auto it = find(types.begin(), types.end(), target);
        if (it != types.end()) {
            idx.p2 = distance(types.begin(), it);
            idx.c2 = target;
            break;
        }
    }
    // Fallback P2
    if (idx.p2 == -1) {
        for (size_t i=0; i<types.size(); i++) {
            if((int)i == idx.p1) continue; // 不重复
            if(types[i][0]!='C') continue;
            if((sys=='G'||sys=='R') && types[i][1]=='2') { idx.p2=i; idx.c2=types[i]; break; }
            if(sys=='E' && (types[i][1]=='5'||types[i][1]=='7')) { idx.p2=i; idx.c2=types[i]; break; }
            if(sys=='C' && (types[i][1]=='6'||types[i][1]=='7')) { idx.p2=i; idx.c2=types[i]; break; }
        }
    }
    // 查找 D2 (对应 P2)
    if (idx.p2 != -1) {
        string tD = idx.c2; tD[0] = 'D';
        auto it = find(types.begin(), types.end(), tD);
        if (it != types.end()) idx.d2 = distance(types.begin(), it);
    }
    return idx;
}

vector<vector<ObsData>> readObsFile(const string& filename, double* approxPos) {
    vector<vector<ObsData>> allEpochs;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open obs file " << filename << endl;
        return allEpochs;
    }

    string line;
    map<char, ColIdx> sys_col_idx;

    // --- Header 处理 ---
    while (getline(file, line)) {
        if (line.length() < 60) continue;

        // 读取概略坐标
        if (line.find("APPROX POSITION XYZ") != string::npos) {
            approxPos[0] = str2num(line.substr(0, 14));
            approxPos[1] = str2num(line.substr(14, 14));
            approxPos[2] = str2num(line.substr(28, 14));
        }

        // 读取观测类型
        if (line.find("SYS / # / OBS TYPES") != string::npos) {
            char sys = line[0];
            if (sys == ' ') continue; // 有时第一列为空 (Legacy RINEX)
            try {
                int n_obs = stoi(line.substr(3, 3));
                vector<string> types;
                // 读取当前行 (最多13个)
                for(int i=0; i<13 && i<n_obs; i++) {
                    if (line.length() >= 7 + i*4 + 3)
                        types.push_back(line.substr(7 + i*4, 3));
                }
                // 如果超过13个，需读续行
                int lines_to_read = (n_obs - 1) / 13;
                for(int k=0; k<lines_to_read; k++) {
                    getline(file, line);
                    for(int i=0; i<13; i++) {
                        int idx = (k+1)*13 + i;
                        if (idx >= n_obs) break;
                        if (line.length() >= 7 + i*4 + 3)
                            types.push_back(line.substr(7 + i*4, 3));
                    }
                }
                sys_col_idx[sys] = map_signal_col(sys, types);
            } catch (...) {}
        }

        if (line.find("END OF HEADER") != string::npos) break;
    }
    cout << "------------------------------------------" << endl;

    // --- Body 处理 ---
    vector<ObsData> epochObs;
    GPSTime epochTime = {0, 0};

    while (getline(file, line)) {
        if (line.empty()) continue;

        // 历元头 >
        if (line[0] == '>') {
            // 保存上一历元
            if (!epochObs.empty()) {
                allEpochs.push_back(epochObs);
                epochObs.clear();
            }
            try {
                int yr = stoi(line.substr(2, 4));
                int mon = stoi(line.substr(7, 2));
                int day = stoi(line.substr(10, 2));
                int hr = stoi(line.substr(13, 2));
                int min = stoi(line.substr(16, 2));
                double sec = str2num(line.substr(19, 11));

                double ep[6] = {(double)yr, (double)mon, (double)day, (double)hr, (double)min, sec};
                epochTime = epoch2time(ep);
            } catch(...) { continue; }
        } else {
            // 观测数据行
            char sys = line[0];
            if (sys_col_idx.find(sys) == sys_col_idx.end()) continue;

            try {
                ObsData obs = {0};
                obs.time = epochTime;
                obs.sys = sys;
                obs.sat = stoi(line.substr(1, 2));

                ColIdx idx = sys_col_idx[sys];

                // Lambda helper to read val
                auto read_val = [&](int col) -> double {
                    if (col < 0) return 0.0;
                    int start = 3 + col * 16;
                    if (line.length() < start + 14) return 0.0; // 简单边界检查
                    return str2num(line.substr(start, 14));
                };

                obs.P1 = read_val(idx.p1);
                obs.P2 = read_val(idx.p2);
                obs.D1 = read_val(idx.d1);
                obs.D2 = read_val(idx.d2);

                obs.codeP1 = idx.c1;
                obs.codeP2 = idx.c2;

                // 只有当有数据时才添加
                if (obs.P1 != 0.0 || obs.D1 != 0.0) {
                    epochObs.push_back(obs);
                }
            } catch (...) {}
        }
    }
    // 添加最后一历元
    if (!epochObs.empty()) allEpochs.push_back(epochObs);

    cout << "Loaded " << allEpochs.size() << " epochs." << endl;
    return allEpochs;
}