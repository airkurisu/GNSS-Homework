/*
 * src/gtime.cpp - 基于 TimeConverter 逻辑的时间处理实现
 */

#include "gtime.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdlib>

using namespace std;

const double LEAPS = 18.0; // 当前跳秒

// --- 移植自你提供的代码: 基础时间辅助函数 ---

// 判断闰年
static bool isLeapYear(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

// 获取该日期在当年是第几天 (Day of Year)
static int getDayOfYear(int year, int month, int day) {
    int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (isLeapYear(year)) {
        days_in_month[1] = 29;
    }

    int doy = 0;
    for (int i = 0; i < month - 1; i++) {
        doy += days_in_month[i];
    }
    doy += day;
    return doy;
}

// 计算两个日期之间的天数差 (核心逻辑)
// 移植自 TimeConverter::daysBetween
static int daysBetween(int year1, int month1, int day1, int year2, int month2, int day2) {
    // 简单计算：年份 * 365 + 该年的第几天
    // 注意：这里用 long long 防止溢出，虽然 int 在 gnss 范围内通常够用
    long long days1 = (long long)year1 * 365 + getDayOfYear(year1, month1, day1);
    long long days2 = (long long)year2 * 365 + getDayOfYear(year2, month2, day2);

    // 补回闰年多出来的天数
    // 遍历两个年份之间的所有年份
    int y_start = std::min(year1, year2);
    int y_end   = std::max(year1, year2);

    for (int y = y_start; y < y_end; y++) {
        if (isLeapYear(y)) {
            if (y < year2) days2++;
            else days1++; // 这种情况一般不会发生，因为我们是用 1980 算到 2024
        }
    }

    return (int)(days2 - days1);
}

// ----------------------------------------------------------------

// 将日历时间 (年,月,日,时,分,秒) 转为 GPS时间
// 逻辑对应 TimeConverter::convertTime
GPSTime epoch2time(const double *ep) {
    GPSTime t = {0, 0.0};

    int year = (int)ep[0];
    int mon  = (int)ep[1];
    int day  = (int)ep[2];
    int hour = (int)ep[3];
    int min  = (int)ep[4];
    double sec = ep[5];

    // GPS 起始时间: 1980年1月6日
    const int GPS_START_YEAR = 1980;
    const int GPS_START_MONTH = 1;
    const int GPS_START_DAY = 6;

    // 1. 计算距离 GPS 起始日期的天数
    int days_from_gps_start = daysBetween(GPS_START_YEAR, GPS_START_MONTH, GPS_START_DAY,
                                          year, mon, day);

    // 2. 计算 GPS 周
    t.week = days_from_gps_start / 7;

    // 3. 计算周内秒 (由天数余数 + 当天时分秒组成)
    // 注意：days_from_gps_start % 7 得到的是这周过了几天 (0~6)
    int day_of_week = days_from_gps_start % 7;

    // 处理负数日期的情况 (虽然一般处理2024年数据不会遇到)
    if (day_of_week < 0) {
        day_of_week += 7;
        t.week -= 1;
    }

    t.sec = day_of_week * 86400.0 + hour * 3600.0 + min * 60.0 + sec;

    return t;
}

// 将 GPS时间 转为 日历时间 (保持原有的高精度反算逻辑，因为 daysBetween 不好逆向)
// 使用标准的 MJD 算法进行反算，这是数学上严谨的
void time2epoch(GPSTime t, double *ep) {
    // 1. 转为儒略日 (JD)
    // GPS起始历元 (1980-01-06) 的 JD = 2444244.5
    double jd = t.week * 7.0 + t.sec / 86400.0 + 2444244.5;

    // 2. JD -> YMD (经典天文算法)
    double z = floor(jd + 0.5);
    double f = jd + 0.5 - z;

    double alpha = floor((z - 1867216.25) / 36524.25);
    double a = z + 1 + alpha - floor(alpha / 4.0);
    double b = a + 1524;
    double c = floor((b - 122.1) / 365.25);
    double d = floor(365.25 * c);
    double e = floor((b - d) / 30.6001);

    ep[2] = b - d - floor(30.6001 * e) + f; // Day
    ep[1] = e < 14 ? e - 1 : e - 13;        // Month
    ep[0] = ep[1] > 2 ? c - 4716 : c - 4715;// Year

    // 提取时分秒
    double day_int = floor(ep[2]);
    double day_frac = ep[2] - day_int;
    ep[2] = day_int;

    double total_sec = day_frac * 86400.0;
    ep[3] = floor(total_sec / 3600.0);       // Hour
    total_sec -= ep[3] * 3600.0;
    ep[4] = floor(total_sec / 60.0);         // Min
    ep[5] = total_sec - ep[4] * 60.0;        // Sec

    // 精度保护
    if (ep[5] >= 60.0) { ep[5] = 0.0; ep[4] += 1; }
    if (ep[4] >= 60.0) { ep[4] = 0.0; ep[3] += 1; }
}

double timediff(GPSTime t1, GPSTime t2) {
    return (t1.week - t2.week) * 604800.0 + (t1.sec - t2.sec);
}

GPSTime gpst2utc(GPSTime t) {
    GPSTime tu = t;
    tu.sec -= LEAPS;
    if (tu.sec < 0) { tu.sec += 604800; tu.week--; }
    return tu;
}

GPSTime utc2gpst(GPSTime t) {
    GPSTime tg = t;
    tg.sec += LEAPS;
    if (tg.sec >= 604800) { tg.sec -= 604800; tg.week++; }
    return tg;
}

double str2num(const char *s, int i, int n) {
    char str[256], *p = str;
    if (i < 0 || n < 1 || sizeof(str) < n + 1) return 0.0;
    for (s += i; *s && --n >= 0; s++) *p++ = *s == 'D' || *s == 'd' ? 'E' : *s;
    *p = '\0';
    return atof(str);
}