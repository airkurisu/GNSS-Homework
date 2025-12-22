#ifndef MYGNSS10_GTIME_H
#define MYGNSS10_GTIME_H

#include <time.h>

struct GPSTime {
    int week;    // GPS周
    double sec;  // 周内秒
};

GPSTime epoch2time(const double *ep);
void time2epoch(GPSTime t, double *ep);
double timediff(GPSTime t1, GPSTime t2);
GPSTime gpst2utc(GPSTime t);
GPSTime utc2gpst(GPSTime t);
double str2num(const char *s, int i, int n);

#endif // GNSS_SPP_GTIME_H