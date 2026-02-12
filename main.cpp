/*
* main.cpp - 项目入口文件
 */
#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include "rinex.h"
#include "spp.h"
#include "spp_gnss.h"

using namespace std;

int main() {
    int mode = 2;
    //if (mode==1)
    //    freopen("out_gps.txt", "w", stdout);
    //else
    //    freopen("out_multi-gnss.txt", "w", stdout);

    // 简单交互菜单

    cout << "Select Processing Mode:" << endl;
    cout << "1. GPS-Only SPP " << endl;
    cout << "2. Multi-GNSS SPP (GPS/BDS/GAL)" << endl;
    cout << "Enter choice [1-2]: ";

    string nav_file = "data/brdc1590.24p";
    string obs_file = "data/jfng1590.24o";


    cout << "Step 1: Reading Navigation file..." << endl;
    auto navList = readNavFile(nav_file);
    if (navList.empty()) return 1;

    cout << "Step 2: Reading Observation file..." << endl;
    double refPos[3] = {0};
    auto obsList = readObsFile(obs_file, refPos);
    if (obsList.empty()) return 1;

    // refPos[0] = ...; // 若需手动指定参考坐标

    // 分支调用
    if (mode == 2) {
        spp_gnss_process(obsList, navList, refPos);
    } else {
        spp_process(obsList, navList, refPos);
    }

    return 0;
}