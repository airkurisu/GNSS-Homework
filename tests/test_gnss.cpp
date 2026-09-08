#include "gtime.h"
#include "matrix.h"
#include "orbit.h"
#include "rinex.h"
#include "spp.h"
#include "tropo.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
constexpr double PI=3.14159265358979323846;
int checks=0;
void check(bool value, const char* message) {
    ++checks;
    if (!value) throw std::runtime_error(message);
}
void near(double value, double expected, double tolerance, const char* message) {
    if (!std::isfinite(value) || std::fabs(value-expected)>tolerance) {
        std::cerr<<message<<": actual="<<std::setprecision(15)<<value<<", expected="<<expected<<'\n';
        check(false,message);
    }
    ++checks;
}
template<class F> void rejects(F action, const char* message) {
    bool threw=false;
    try { action(); } catch (const std::exception&) { threw=true; }
    check(threw,message);
}
std::string header(std::string body, const std::string& label) {
    body.resize(60,' '); return body+label+"\n";
}
std::string obsField(double value) {
    std::ostringstream s; s<<std::fixed<<std::setprecision(3)<<std::setw(14)<<value<<"  "; return s.str();
}
std::string navField(double value) {
    std::ostringstream s; s<<std::scientific<<std::setprecision(12)<<std::setw(19)<<value; return s.str();
}
void write(const std::string& path, const std::string& text) { std::ofstream f(path); f<<text; }
void timeTests() {
    double ep[]={1980,1,6,0,0,0};
    auto t=epoch2time(ep); check(t.week==0 && t.sec==0,"GPS origin");
    for (int year:{1979,1980,1999,2000,2024,2100}) {
        double date[]={double(year),3,1,23,59,59.1234567}, back[6];
        auto a=epoch2time(date); time2epoch(a,back);
        for (int i=0; i<6; ++i) near(back[i],date[i],1e-9,"Calendar round trip");
    }
    auto crossed=timeadd({2300,0.1},-0.2);
    check(crossed.week==2299,"Negative week rollover"); near(crossed.sec,604799.9,1e-9,"Rollover seconds");
    double a[]={2016,1,1,0,0,0}, b[]={2024,1,1,0,0,0};
    near(timediff(utc2gpst(epoch2time(a)),epoch2time(a)),17,0,"Historical leap seconds");
    near(timediff(utc2gpst(epoch2time(b)),epoch2time(b)),18,0,"2024 leap seconds");
    near(timediff(gpst2utc(utc2gpst(epoch2time(a))),epoch2time(a)),0,0,"UTC round trip");
    double invalid[]={2023,2,29,0,0,0}; rejects([&]{epoch2time(invalid);},"Reject invalid date");
    near(str2num("1",50,19),0,0,"Short numeric string must not overread");
    near(str2num(" 1.25D-03",0,10),0.00125,1e-15,"Fortran exponent");
}
void geometryTests() {
    const double f=1/298.257223563, e2=f*(2-f), lat=30*PI/180, lon=115*PI/180, height=123.456;
    const double n=6378137/std::sqrt(1-e2*std::sin(lat)*std::sin(lat));
    double r[]={(n+height)*std::cos(lat)*std::cos(lon),(n+height)*std::cos(lat)*std::sin(lon),
                (n*(1-e2)+height)*std::sin(lat)}, pos[3];
    ecef2pos(r,pos);
    near(pos[0],lat,1e-12,"Geodetic latitude"); near(pos[1],lon,1e-12,"Longitude");
    near(pos[2],height,1e-6,"Ellipsoid height");
    double up[]={std::cos(lat)*std::cos(lon),std::cos(lat)*std::sin(lon),std::sin(lat)}, enu[3];
    ecefVectorToEnu(up,r,enu);
    near(enu[0],0,1e-12,"ENU east"); near(enu[1],0,1e-12,"ENU north"); near(enu[2],1,1e-12,"ENU up");
    double pole[]={0,0,6378137*(1-f)+200}; ecef2pos(pole,pos);
    near(pos[0],PI/2,1e-12,"North pole latitude"); near(pos[2],200,1e-7,"Polar height");
    double origin[3]{}; ecef2pos(origin,pos); check(std::isfinite(pos[2]),"Origin remains finite");
    double H[]={1,0,0,1,1,1,2,1}, v[]={2,3,5,7}, w[]={1,1,1,1}, x[2];
    check(leastSquares(H,v,w,4,2,x),"QR solve"); near(x[0],2,1e-12,"QR x"); near(x[1],3,1e-12,"QR y");
    double singular[]={1,1,2,2,3,3,4,4}; check(!leastSquares(singular,v,w,4,2,x),"Rank deficient QR");
    double nearly[]={1,1,1,1+1e-14}; check(matinv(nearly,2)==-1,"Near singular inverse");
    double m[]={4,7,2,6}; check(matinv(m,2)==0,"Matrix inverse"); near(m[0],0.6,1e-12,"Inverse value");
    double azel[]={0,PI/2}, trop=0,var=0, ground[]={lat,lon,0};
    check(trop_model_prec({2300,0},ground,azel,&trop,&var)==1 && trop>2 && trop<3,"Zenith atmosphere");
}
void orbitTests() {
    NavData n; n.sys='G'; n.sat=1; n.sqrtA=5153.7955; n.e=0.02; n.i0=0.94;
    n.M0=0.3; n.Omega0=1.2; n.omega=0.7; n.OmegaDot=-8e-9;
    n.toe={2300,100000}; n.toeSow=100000; n.toc={2300,99000};
    n.a0=2e-4; n.a1=1e-10; n.a2=2e-15;
    for (char sys:{'G','E','C'}) for (int prn:{1,20}) {
        n.sys=sys; n.sat=prn;
        GPSTime t={2300,100123}; double p[3],v[3],c,d,plus[3],minus[3],cp,cm;
        satPosVel(t,n,p,v,&c,&d);
        satPos(timeadd(t,0.01),n,plus,&cp); satPos(timeadd(t,-0.01),n,minus,&cm);
        for (int j=0; j<3; ++j) near(v[j],(plus[j]-minus[j])/0.02,2e-5,"Analytic orbit velocity");
        near(d,(cp-cm)/0.02,1e-17,"Clock drift includes relativity derivative");
    }
    n.sys='G'; n.e=0;
    double p[3],c; satPos({2300,100123},n,p,&c);
    near(c,n.a0+n.a1*1123+n.a2*1123*1123,1e-16,"Clock polynomial must use TOC");
}
void parserTests() {
    const std::string path="test-observation.rnx", navPath="test-navigation.rnx";
    std::string base=header("     3.04           O","RINEX VERSION / TYPE")+
        header("G    8 C1W C1C C2W C2X D1W D1C D2W D2X","SYS / # / OBS TYPES")+
        header("","END OF HEADER");
    // Preferred codes missing on this satellite: choose alternatives per observation.
    std::string data="G01"+std::string(16,' ')+obsField(22000000)+std::string(16,' ')+obsField(22000005)+
        std::string(16,' ')+obsField(0)+std::string(16,' ')+obsField(0)+"\n";
    write(path,base+"> 2024 06 07 00 00 00.0000000  0  1\n"+data+
        "> 2024 06 07 00 00 01.0000000  5  1\nexternal event\n"+
        "> 2024 06 07 00 00 30.0000000  0  1\nG01"+obsField(22000000)+"\n");
    double reference[3]; auto o=readObsFile(path,reference);
    check(o.size()==2,"Event records are not observations");
    check(o[0][0].codeP1=="C1C" && o[0][0].codeP2=="C2X","Per-satellite signal fallback");
    check(o[0][0].D1==0 && o[0][0].D2==0,"Measured zero Doppler is valid");
    check(o[1][0].sat==0 && std::isnan(o[1][0].D1),"Missing data keeps epoch without false zero Doppler");
    write(path,base+"> 2024 06 07 00 00 00.0000000  0  2\n"+data);
    rejects([&]{readObsFile(path,reference);},"Reject truncated epoch");
    write(path,base+"> 2024 06 07 00 00 00.0000000  0  1\nG01           nan\n");
    rejects([&]{readObsFile(path,reference);},"Reject nonfinite observation");
    write(path,base+"> 2024 06 07 00 00 00.0000000  4  0\n");
    rejects([&]{readObsFile(path,reference);},"Explicit unsupported header update");
    std::string nb=header("     3.04           N","RINEX VERSION / TYPE")+header("","END OF HEADER");
    std::string record="C01 2024 06 08 23 59 50"+navField(0)+navField(0)+navField(0)+"\n";
    double rows[7][4]={{0,0,0,0},{0,0.01,0,5153.7955},{604790,0,0,0},{0.94,0,0,0},
                       {0,0,961,0},{0,0,1e-8,2e-8},{0,0,0,0}};
    for (auto& row:rows) { record+="    "; for (double x:row) record+=navField(x); record+='\n'; }
    write(navPath,nb+record); auto nav=readNavFile(navPath);
    check(nav.size()==1 && nav[0].toe.week==2318,"BDT to GPST week rollover");
    near(nav[0].toe.sec,4,1e-9,"BDT to GPST seconds");
    near(timediff(nav[0].toe,nav[0].toc),0,1e-9,"BDT TOC and TOE consistency");
    near(nav[0].toeSow,604790,0,"Native TOE seconds preserved");
    near(nav[0].tgd2,2e-8,1e-18,"Second BDS group delay");
    write(navPath,nb+record.substr(0,record.find('\n')+1));
    rejects([&]{readNavFile(navPath);},"Reject truncated navigation");
    std::remove(path.c_str()); std::remove(navPath.c_str());
}

void syntheticTests() {
    // Known moving receiver and large, independent system clock drifts. This catches
    // clock-sign errors that look harmless on a nearly stationary real receiver.
    const double rr[]={6378137,0,0}, vr[]={12,-3,1.5};
    GPSTime time={2300,100000};
    std::vector<NavData> nav;
    std::vector<ObsData> obs;
    for (char sys:{'G','C'}) {
        int count=0;
        for (int k=0; k<120 && count<8; ++k) {
            NavData n; n.sys=sys; n.sat=10+count; n.sqrtA=5153.7955; n.e=0.01;
            n.i0=0.96; n.Omega0=(k%12)*2*PI/12; n.M0=(k/12)*2*PI/10;
            n.toe=time; n.toeSow=sys=='C'?time.sec-14:time.sec; n.toc=timeadd(time,-1000);
            n.a0=1e-4; n.a1=(count+1)*1e-10; n.a2=1e-15;
            n.tgd=sys=='C'?1e-8:0;
            double rs[3],vs[3],c,d;
            satPosVel(time,n,rs,vs,&c,&d);
            if (rs[0]-rr[0]<0.3*std::hypot(std::hypot(rs[0]-rr[0],rs[1]),rs[2])) continue;
            ObsData o; o.sys=sys; o.sat=n.sat; o.time=time;
            o.codeP1=sys=='G'?"C1W":"C2I"; o.codeP2=sys=='G'?"C2W":"C6I";
            double f1=0,f2=0; check(observationFrequencies(o,n,f1,f2),"Synthetic frequencies");
            const double bias=sys=='G'?30000:-20000, drift=sys=='G'?20:-30;
            o.P1=o.P2=24000000;
            for (int iter=0; iter<8; ++iter) {
                GPSTime tx=timeadd(time,-o.P1/CLIGHT);
                satPos(tx,n,rs,&c); tx=timeadd(tx,-c);
                satPosVel(tx,n,rs,vs,&c,&d);
                // Independent exact rotation instead of the solver's first-order range formula.
                double rho=std::hypot(std::hypot(rs[0]-rr[0],rs[1]),rs[2]);
                double theta=OMEGA_EARTH*rho/CLIGHT;
                double rotatedX=std::cos(theta)*rs[0]+std::sin(theta)*rs[1];
                double rotatedY=-std::sin(theta)*rs[0]+std::cos(theta)*rs[1];
                rho=std::hypot(std::hypot(rotatedX-rr[0],rotatedY),rs[2]);
                double azel[]={0,std::atan2(rs[0]-rr[0],std::hypot(rs[1],rs[2]))};
                double pos[]={0,0,0}, trop=0,var=0; trop_model_prec(time,pos,azel,&trop,&var);
                double p=rho+bias-CLIGHT*c+trop;
                o.P1=p+4+CLIGHT*n.tgd;
                o.P2=p+4*f1*f1/(f2*f2);
            }
            double line[]={rs[0]-rr[0],rs[1],rs[2]}, length=std::hypot(std::hypot(line[0],line[1]),line[2]);
            double rate=0; for (int j=0; j<3; ++j) rate+=(vs[j]-vr[j])*line[j]/length;
            rate+=OMEGA_EARTH/CLIGHT*(vs[1]*rr[0]+rs[1]*vr[0]-vs[0]*rr[1]-rs[0]*vr[1]);
            rate+=drift-CLIGHT*d;
            o.D1=-rate*f1/CLIGHT; o.D2=-rate*f2/CLIGHT;
            nav.push_back(n); obs.push_back(o); ++count;
        }
        check(count==8,"Enough synthetic satellites");
    }
    SPPSolver solver(nav); auto r=solver.solve(obs,rr);
    check(r.positionValid && r.velocityValid,"Synthetic combined solution");
    for (int j=0; j<3; ++j) {
        near(r.pos[j],rr[j],0.002,"Synthetic receiver position");
        near(r.vel[j],vr[j],2e-5,"Synthetic receiver velocity");
    }
    near(r.clockDrift.at('G')*CLIGHT,20,2e-5,"GPS clock drift sign");
    near(r.clockDrift.at('C')*CLIGHT,-30,2e-5,"BDS clock drift sign");
    auto gross=obs;
    gross[0].P1+=150; gross[0].P2+=150;
    auto screened=solver.solve(gross,rr);
    check(screened.positionValid && screened.nRejected==1,"Remove gross pseudorange error using SPP residuals");
    for (int j=0; j<3; ++j) near(screened.pos[j],rr[j],0.003,"Position after outlier screening");
    auto zeroDoppler=obs;
    for (auto& o:zeroDoppler) o.D1=o.D2=0;
    check(solver.solve(zeroDoppler,rr).velocityValid,"Zero Doppler rows participate in velocity solve");
    auto few=obs; few.resize(3);
    check(!solver.solve(few,rr).positionValid,"Three satellites cannot determine a position");
    double zero[3]{}; auto cold=solver.solve(obs,zero);
    check(cold.positionValid,"Synthetic cold start");
    for (int j=0; j<3; ++j) near(cold.pos[j],rr[j],0.002,"Cold-start position");
    // An observed system with no usable navigation must not create an empty clock column.
    ObsData missing=obs[0]; missing.sys='E'; missing.sat=1; obs.push_back(missing);
    check(solver.solve(obs,rr).velocityValid,"Absent system excluded from design matrix");
    for (auto& o:obs) o.D1=o.D2=NAN;
    auto noVelocity=solver.solve(obs,rr);
    check(noVelocity.positionValid && !noVelocity.velocityValid,"Position retained without Doppler");
    auto stale=nav; for (auto& n:stale) n.toe=timeadd(n.toe,-86400);
    check(!SPPSolver(stale).solve(obs,rr).positionValid,"Reject stale ephemerides");
    for (auto& n:nav) n.health=1;
    check(!SPPSolver(nav).solve(obs,rr).positionValid,"Reject unhealthy satellites");
}
void realDataTests() {
    auto nav=readNavFile(std::string(GNSS_SOURCE_DIR)+"/data/brdc1590.24p");
    double reference[3]; auto obs=readObsFile(std::string(GNSS_SOURCE_DIR)+"/data/jfng1590.24o",reference);
    check(obs.size()==2880,"Full-day sample epoch count");
    SPPSolver solver(nav);
    for (size_t i=0; i<obs.size(); i+=240) for (const std::string sys:{"G","GCER"}) {
        auto r=solver.solve(obs[i],reference,sys);
        check(r.positionValid && r.velocityValid,"Real sample solution");
        double err=0,speed=0;
        for (int j=0; j<3; ++j) { err+=std::pow(r.pos[j]-reference[j],2); speed+=r.vel[j]*r.vel[j]; }
        check(std::sqrt(err)<20,"Real position sanity bound, not a truth accuracy claim");
        check(std::sqrt(speed)<0.5,"Static receiver velocity sanity bound");
    }
}
}
int main() {
    try {
        timeTests(); geometryTests(); orbitTests(); parserTests(); syntheticTests(); realDataTests();
        std::cout<<"Passed "<<checks<<" checks\n"; return 0;
    } catch (const std::exception& e) { std::cerr<<"FAILED: "<<e.what()<<'\n'; return 1; }
}
