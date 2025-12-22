/*
 * main.cpp - 项目入口文件
 * 项目: MYGNSS10
 * 功能: 读取数据并调用SPP处理模块
 */
#include <iostream>
#include <vector>
#include<iomanip>
#include <string>
#include "rinex.h"
#include "spp.h"

using namespace std;

int main() {
    freopen("out.txt", "w", stdout);
    string nav_file = "data/brdc1590.24p";
    string obs_file = "data/jfng1590.24o";

    cout << "==================================================" << endl;
    cout << "       MYGNSS10: SPP & Velocity Estimator         " << endl;
    cout << "==================================================" << endl;

    // --- 2. 读取广播星历 ---
    cout << "Step 1: Reading Navigation file..." << endl;
    cout << "Path: " << nav_file << endl;

    auto navList = readNavFile(nav_file);
    if (navList.empty()) {
        cerr << "[Error] Failed to read navigation data or file is empty." << endl;
        return 1;
    }
    cout << "Done. Loaded " << navList.size() << " ephemeris records." << endl << endl;

    // --- 3. 读取观测文件 ---
    cout << "Step 2: Reading Observation file..." << endl;
    cout << "Path: " << obs_file << endl;

    // 用于存储从 RINEX 头文件中读取的近似坐标 (XYZ)
    // 这通常作为参考真值 (Reference Position) 用于计算 ENU 误差
    double refPos[3] = {0};

    auto obsList = readObsFile(obs_file, refPos);
    if (obsList.empty()) {
        cerr << "[Error] Failed to read observation data or file is empty." << endl;
        return 1;
    }
    cout << "Done. Loaded " << obsList.size() << " epochs." << endl;
    cout << "Reference Position (from Header): "
         << fixed << setprecision(3)
         << refPos[0] << ", " << refPos[1] << ", " << refPos[2] << endl << endl;

    // --- 4. 执行 SPP 定位与测速 ---
    cout << "Step 3: Processing..." << endl;

    // 如果观测文件头中没有坐标，您可以手动在此处指定真值 (例如 IGS 站精确坐标)
    // refPos[0] = -2267749.0; refPos[1] = 5009154.0; refPos[2] = 3221290.0;

    // 调用核心处理函数
    spp_process(obsList, navList, refPos);

    cout << endl << "Program finished." << endl;
    return 0;
}