#include "gtime.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace {
bool leap(int y) { return y % 4 == 0 && (y % 100 != 0 || y % 400 == 0); }
int monthDays(int y, int m) {
    const int days[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    return days[m-1] + (m == 2 && leap(y));
}
// Published UTC leap-second effective dates, through 2017-01-01.
const int leaps[][3] = {
    {1981,7,1},{1982,7,1},{1983,7,1},{1985,7,1},{1988,1,1},{1990,1,1},
    {1991,1,1},{1992,7,1},{1993,7,1},{1994,7,1},{1996,1,1},{1997,7,1},
    {1999,1,1},{2006,1,1},{2009,1,1},{2012,7,1},{2015,7,1},{2017,1,1}
};
}
GPSTime timeadd(GPSTime t, double seconds) {
    t.sec += seconds;
    const int weeks = static_cast<int>(std::floor(t.sec / 604800.0));
    t.week += weeks;
    t.sec -= weeks * 604800.0;
    return t;
}
GPSTime epoch2time(const double* ep) {
    for (int i=0; i<6; ++i) if (!std::isfinite(ep[i])) throw std::runtime_error("Non-finite epoch");
    for (int i=0; i<5; ++i) if (ep[i]!=std::floor(ep[i])) throw std::runtime_error("Non-integral calendar field");
    if (ep[0]<1900 || ep[0]>2199 || ep[1]<1 || ep[1]>12) throw std::runtime_error("Invalid year/month");
    if (ep[2]<1 || ep[2]>31) throw std::runtime_error("Invalid day");
    const int y=static_cast<int>(ep[0]), m=static_cast<int>(ep[1]), d=static_cast<int>(ep[2]);
    if (d<1 || d>monthDays(y,m) || ep[3]<0 || ep[3]>=24 || ep[4]<0 || ep[4]>=60 || ep[5]<0 || ep[5]>=60)
        throw std::runtime_error("Invalid calendar epoch (leap-second labels unsupported)");
    int days=-5;
    for (int k=1980; k<y; ++k) days += leap(k) ? 366 : 365;
    for (int k=y; k<1980; ++k) days -= leap(k) ? 366 : 365;
    for (int k=1; k<m; ++k) days += monthDays(y,k);
    days += d-1;
    const int week=static_cast<int>(std::floor(days/7.0));
    return timeadd({week,(days-week*7)*86400.0+ep[3]*3600+ep[4]*60+ep[5]},0);
}
void time2epoch(GPSTime t, double* ep) {
    t=timeadd(t,0);
    int days=t.week*7+static_cast<int>(std::floor(t.sec/86400))+5;
    double sec=t.sec-std::floor(t.sec/86400)*86400;
    int y=1980, m=1;
    while (days<0) { --y; days+=leap(y)?366:365; }
    while (days>=(leap(y)?366:365)) { days-=leap(y)?366:365; ++y; }
    while (days>=monthDays(y,m)) { days-=monthDays(y,m); ++m; }
    ep[0]=y; ep[1]=m; ep[2]=days+1;
    ep[3]=std::floor(sec/3600); sec-=ep[3]*3600;
    ep[4]=std::floor(sec/60); ep[5]=sec-ep[4]*60;
}
double timediff(GPSTime a, GPSTime b) { return (a.week-b.week)*604800.0+a.sec-b.sec; }
GPSTime utc2gpst(GPSTime t) {
    int offset=0;
    for (const auto& date:leaps) {
        double ep[]={double(date[0]),double(date[1]),double(date[2]),0,0,0};
        if (timediff(t,epoch2time(ep))>=0) ++offset;
    }
    return timeadd(t,offset);
}
GPSTime gpst2utc(GPSTime t) {
    for (int i=17; i>=0; --i) {
        double ep[]={double(leaps[i][0]),double(leaps[i][1]),double(leaps[i][2]),0,0,0};
        GPSTime utc=timeadd(t,-(i+1));
        if (timediff(utc,epoch2time(ep))>=0) return utc;
    }
    return timeadd(t,0);
}
double str2num(const char* s, int i, int n) {
    if (!s || i<0 || n<1 || static_cast<size_t>(i)>=std::strlen(s)) return 0;
    std::string field=std::string(s).substr(static_cast<size_t>(i),static_cast<size_t>(n));
    for (char& c:field) if (c=='D' || c=='d') c='E';
    return std::strtod(field.c_str(),nullptr);
}
