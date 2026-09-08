#include "spp.h"
#include "matrix.h"
#include "orbit.h"
#include "tropo.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>

namespace {
constexpr double PI=3.14159265358979323846;
constexpr double A_WGS=6378137.0, F_WGS=1.0/298.257223563;
constexpr double MASK=10*PI/180;
double norm3(const double* a) { return std::hypot(std::hypot(a[0],a[1]),a[2]); }
bool finite3(const double* a) { return std::isfinite(a[0]) && std::isfinite(a[1]) && std::isfinite(a[2]); }
double ionoFree(double a, double b, double f1, double f2) {
    return (f1*f1*a-f2*f2*b)/(f1*f1-f2*f2);
}
bool validPosition(const double* r) { return finite3(r) && norm3(r)>A_WGS/2; }
struct Measurement {
    ObsData obs;
    double rs[3]{}, vs[3]{}, clock=0, drift=0, pseudorange=0, rangeRate=NAN;
    double codeSigma=1; // IF noise amplification for an assumed 1 m noise per raw code
};
struct Row {
    char sys;
    int sat=0;
    double h[3];
    double residual, weight;
    double sigma=1;
};
std::vector<char> systemsIn(const std::vector<Row>& rows) {
    std::set<char> systems;
    for (const auto& row:rows) systems.insert(row.sys);
    return {systems.begin(),systems.end()};
}
bool fit(const std::vector<Row>& rows, const std::vector<char>& systems,
         std::vector<double>& dx, double* pdop=nullptr) {
    const int n=static_cast<int>(rows.size()), m=3+static_cast<int>(systems.size());
    if (n<m) return false;
    std::vector<double> H, v, w;
    for (const auto& row:rows) {
        H.insert(H.end(),row.h,row.h+3);
        for (char sys:systems) H.push_back(sys==row.sys?1:0);
        v.push_back(row.residual); w.push_back(row.weight);
    }
    dx.assign(m,0);
    if (!leastSquares(H.data(),v.data(),w.data(),n,m,dx.data())) return false;
    if (pdop) {
        std::vector<double> q(m*m,0);
        for (int i=0; i<n; ++i) for (int j=0; j<m; ++j) for (int k=0; k<m; ++k)
            q[j*m+k]+=H[i*m+j]*H[i*m+k];
        if (matinv(q.data(),m)!=0) return false;
        *pdop=std::sqrt(q[0]+q[m+1]+q[2*m+2]);
        if (!std::isfinite(*pdop)) return false;
    }
    return true;
}
int outlier(const std::vector<Row>& rows, const std::vector<char>& systems) {
    const int m=3+static_cast<int>(systems.size());
    // A minimally determined solution cannot support an independent residual check.
    if (static_cast<int>(rows.size())<m+2) return -1;
    std::vector<double> H, q(m*m,0);
    for (const auto& row:rows) {
        H.insert(H.end(),row.h,row.h+3);
        for (char sys:systems) H.push_back(sys==row.sys?1:0);
    }
    for (size_t i=0; i<rows.size(); ++i) for (int j=0; j<m; ++j) for (int k=0; k<m; ++k)
        q[j*m+k]+=H[i*m+j]*rows[i].weight*H[i*m+k];
    if (matinv(q.data(),m)!=0) return -1;
    double worst=6.0; // deliberately conservative six-sigma gross-error screen
    int index=-1;
    for (size_t i=0; i<rows.size(); ++i) {
        // Row noise differs by IF frequency pair. Compute the actual residual
        // variance using the linear weighted projection, not just (1-leverage).
        double variance=0;
        for (size_t k=0; k<rows.size(); ++k) {
            double projection=0;
            for (int a=0; a<m; ++a) for (int b=0; b<m; ++b)
                projection+=H[i*m+a]*q[a*m+b]*H[k*m+b]*rows[k].weight;
            const double coefficient=(i==k?1.0:0.0)-projection;
            variance+=coefficient*coefficient*rows[k].sigma*rows[k].sigma/rows[k].weight;
        }
        if (variance<1e-8) continue; // a singleton system clock absorbs its own row
        double score=std::fabs(rows[i].residual)/std::sqrt(variance);
        if (score>worst) { worst=score; index=static_cast<int>(i); }
    }
    return index;
}
double elevation(const double* rr, const double* rs) {
    double delta[3],enu[3];
    for (int i=0; i<3; ++i) delta[i]=rs[i]-rr[i];
    ecefVectorToEnu(delta,rr,enu);
    return std::atan2(enu[2],std::hypot(enu[0],enu[1]));
}
double weight(const Measurement& s, double el) {
    double w=std::sin(el)*std::sin(el);
    if (s.obs.sys=='C' && (s.obs.sat<=5 || s.obs.sat>=59)) w*=0.1;
    return w;
}
}

void ecef2pos(const double* r, double* pos) {
    const double e2=F_WGS*(2-F_WGS), p=std::hypot(r[0],r[1]);
    pos[1]=p>1e-12 ? std::atan2(r[1],r[0]) : 0;
    if (p<1e-12) {
        pos[0]=r[2]==0 ? 0 : std::copysign(PI/2,r[2]);
        pos[2]=std::fabs(r[2])-A_WGS*(1-F_WGS);
        return;
    }
    double z=r[2], radius=A_WGS;
    for (int i=0; i<20; ++i) {
        double sinLat=z/std::hypot(p,z);
        radius=A_WGS/std::sqrt(1-e2*sinLat*sinLat);
        double next=r[2]+radius*e2*sinLat;
        if (std::fabs(next-z)<1e-8) { z=next; break; }
        z=next;
    }
    pos[0]=std::atan2(z,p);
    pos[2]=std::hypot(p,z)-radius;
}
void ecefVectorToEnu(const double* v, const double* reference, double* enu) {
    double pos[3]; ecef2pos(reference,pos);
    double sl=std::sin(pos[1]), cl=std::cos(pos[1]), sp=std::sin(pos[0]), cp=std::cos(pos[0]);
    enu[0]=-sl*v[0]+cl*v[1];
    enu[1]=-sp*cl*v[0]-sp*sl*v[1]+cp*v[2];
    enu[2]=cp*cl*v[0]+cp*sl*v[1]+sp*v[2];
}
bool observationFrequencies(const ObsData& o, const NavData& nav, double& f1, double& f2) {
    if (o.codeP1.size()!=3 || o.codeP2.size()!=3) return false;
    char a=o.codeP1[1], b=o.codeP2[1];
    if (o.sys=='G' && a=='1' && b=='2') { f1=1.57542e9; f2=1.22760e9; return true; }
    if (o.sys=='C' && a=='2' && (b=='6' || b=='7')) {
        f1=1.561098e9; f2=b=='6'?1.26852e9:1.20714e9; return true;
    }
    if (o.sys=='E' && a=='1' && (b=='5' || b=='7')) {
        f1=1.57542e9; f2=b=='5'?1.17645e9:1.20714e9; return true;
    }
    if (o.sys=='R' && a=='1' && b=='2' && nav.freq_num>=-7 && nav.freq_num<=6) {
        f1=1.602e9+nav.freq_num*0.5625e6; f2=1.246e9+nav.freq_num*0.4375e6; return true;
    }
    return false;
}
SPPSolver::SPPSolver(const std::vector<NavData>& navigation) {
    for (const auto& n:navigation) navigation_[{n.sys,n.sat}].push_back(n);
}
const NavData* SPPSolver::selectEphemeris(const ObsData& obs) const {
    auto found=navigation_.find({obs.sys,obs.sat});
    if (found==navigation_.end()) return nullptr;
    const double maxAge=obs.sys=='R'?1800:(obs.sys=='G'?7200:(obs.sys=='E'?10800:14400));
    const NavData* best=nullptr; double age=maxAge+1;
    for (const auto& n:found->second) {
        if (n.health!=0) continue;
        // Galileo F/NAV clock uses E1/E5a, I/NAV uses E1/E5b (RINEX bits 8/9).
        if (obs.sys=='E') {
            const int clockBit=obs.codeP2.size()==3 && obs.codeP2[1]=='5'?256:512;
            if (!(n.dataSources&clockBit)) continue;
        }
        double dt=std::fabs(timediff(obs.time,n.toe));
        if (dt<=maxAge && dt<age) { age=dt; best=&n; }
    }
    return best;
}
SPPResult SPPSolver::solve(const std::vector<ObsData>& epoch, const double* initial,
                           const std::string& enabled) const {
    SPPResult result;
    if (epoch.empty()) return result;
    result.time=epoch.front().time;
    if (validPosition(initial)) std::copy(initial,initial+3,result.pos);
    std::vector<Measurement> measurements;
    std::set<std::pair<char,int>> seen;
    for (const auto& obs:epoch) {
        if (enabled.find(obs.sys)==std::string::npos || !std::isfinite(obs.P1) || !std::isfinite(obs.P2) ||
            obs.P1<=0 || obs.P2<=0 || std::fabs(timediff(obs.time,result.time))>1e-6) continue;
        if (!seen.insert({obs.sys,obs.sat}).second) continue;
        const NavData* nav=selectEphemeris(obs);
        double f1=0,f2=0;
        if (!nav || !observationFrequencies(obs,*nav,f1,f2)) continue;
        Measurement s; s.obs=obs;
        const double alpha=f1*f1/(f1*f1-f2*f2), beta=f2*f2/(f1*f1-f2*f2);
        s.codeSigma=std::hypot(alpha,beta);
        GPSTime transmit=timeadd(obs.time,-obs.P1/CLIGHT);
        double clock=0;
        satPos(transmit,*nav,s.rs,&clock);
        transmit=timeadd(transmit,-clock); // convert satellite clock time to system time
        satPosVel(transmit,*nav,s.rs,s.vs,&s.clock,&s.drift);
        if (!finite3(s.rs) || !finite3(s.vs) || !std::isfinite(s.clock+s.drift)) continue;
        s.pseudorange=ionoFree(obs.P1,obs.P2,f1,f2);
        if (obs.sys=='C') {
            double delay=ionoFree(nav->tgd,obs.codeP2[1]=='7'?nav->tgd2:0,f1,f2);
            s.pseudorange-=CLIGHT*delay; // BDS broadcast clock is referenced to B3
        }
        if (std::isfinite(obs.D1) && std::isfinite(obs.D2))
            s.rangeRate=ionoFree(-obs.D1*CLIGHT/f1,-obs.D2*CLIGHT/f2,f1,f2);
        measurements.push_back(s);
    }
    std::map<char,double> clocks;
    result.status="position_not_converged";
    for (int iter=0; iter<15; ++iter) {
        std::vector<Row> rows;
        bool located=validPosition(result.pos);
        double pos[3]{}; if (located) ecef2pos(result.pos,pos);
        for (const auto& s:measurements) {
            double delta[3]; for (int j=0; j<3; ++j) delta[j]=s.rs[j]-result.pos[j];
            double distance=norm3(delta);
            if (distance<1) continue;
            double el=located?elevation(result.pos,s.rs):PI/2;
            if (el<MASK) continue;
            double trop=0,var=0,azel[]={0,el};
            if (located && pos[2]>-1000 && pos[2]<20000) trop_model_prec(result.time,pos,azel,&trop,&var);
            // First-order Sagnac range correction, using ECEF at transmission time.
            double rotation=OMEGA_EARTH/CLIGHT*(s.rs[0]*result.pos[1]-s.rs[1]*result.pos[0]);
            Row row{}; row.sys=s.obs.sys; row.sat=s.obs.sat; row.sigma=s.codeSigma;
            for (int j=0; j<3; ++j) row.h[j]=-delta[j]/distance;
            row.h[0]-=OMEGA_EARTH/CLIGHT*s.rs[1]; row.h[1]+=OMEGA_EARTH/CLIGHT*s.rs[0];
            row.residual=s.pseudorange-(distance+rotation+clocks[row.sys]-CLIGHT*s.clock+trop);
            row.weight=weight(s,el); rows.push_back(row);
        }
        auto systems=systemsIn(rows); std::vector<double> dx;
        if (!fit(rows,systems,dx,&result.pdop)) { result.status="insufficient_position_geometry"; return result; }
        for (int j=0; j<3; ++j) result.pos[j]+=dx[j];
        double maxClock=0;
        for (size_t j=0; j<systems.size(); ++j) { clocks[systems[j]]+=dx[3+j]; maxClock=std::max(maxClock,std::fabs(dx[3+j])); }
        if (!finite3(result.pos)) { result.status="nonfinite_position"; return result; }
        if (norm3(dx.data())<1e-4 && maxClock<1e-4) {
            int bad=outlier(rows,systems);
            if (bad>=0) {
                if (result.nRejected>=5) { result.status="excessive_pseudorange_outliers"; return result; }
                const auto rejected=rows[bad];
                measurements.erase(std::remove_if(measurements.begin(),measurements.end(),[&](const Measurement& s) {
                    return s.obs.sys==rejected.sys && s.obs.sat==rejected.sat;
                }),measurements.end());
                ++result.nRejected;
                clocks.clear();
                std::fill(result.pos,result.pos+3,0);
                if (validPosition(initial)) std::copy(initial,initial+3,result.pos);
                iter=-1; // restart after excluding an observation; capped by nRejected
                continue;
            }
            result.positionValid=true; result.nSat=static_cast<int>(rows.size());
            for (char sys:systems) result.clockBias[sys]=clocks[sys]/CLIGHT;
            break;
        }
    }
    if (!result.positionValid) return result;
    std::vector<Row> rows;
    for (const auto& s:measurements) {
        if (!std::isfinite(s.rangeRate)) continue;
        double el=elevation(result.pos,s.rs);
        if (el<MASK) continue;
        double e[3]; for (int j=0; j<3; ++j) e[j]=s.rs[j]-result.pos[j];
        double distance=norm3(e); for (double& x:e) x/=distance;
        double rate=0; for (int j=0; j<3; ++j) rate+=e[j]*s.vs[j];
        // ECEF range rate and Earth rotation correction (RTKLIB resdop convention).
        rate+=OMEGA_EARTH/CLIGHT*(s.vs[1]*result.pos[0]-s.vs[0]*result.pos[1]);
        Row row{}; row.sys=s.obs.sys;
        for (int j=0; j<3; ++j) row.h[j]=-e[j];
        row.h[0]+=OMEGA_EARTH/CLIGHT*s.rs[1]; row.h[1]-=OMEGA_EARTH/CLIGHT*s.rs[0];
        // Pdot = geometric_rate + receiver_clock_drift - c * satellite_clock_drift.
        row.residual=s.rangeRate-(rate-CLIGHT*s.drift);
        row.weight=weight(s,el); rows.push_back(row);
    }
    auto systems=systemsIn(rows); std::vector<double> dx;
    if (!fit(rows,systems,dx)) { result.status="insufficient_velocity_geometry"; return result; }
    std::copy(dx.begin(),dx.begin()+3,result.vel);
    for (size_t j=0; j<systems.size(); ++j) result.clockDrift[systems[j]]=dx[3+j]/CLIGHT;
    result.velocityValid=true; result.nVel=static_cast<int>(rows.size()); result.status="ok";
    return result;
}
int processObservations(const std::vector<std::vector<ObsData>>& observations,
                         const std::vector<NavData>& navigation, const double* reference,
                         const std::string& systems, bool coldStart) {
    SPPSolver solver(navigation);
    bool haveReference=validPosition(reference);
    std::cout<<"week,tow,status,ns,nv,pdop,x,y,z,de,dn,du,ve,vn,vu,rejected\n";
    int np=0,nv=0;
    for (const auto& epoch:observations) {
        double zero[3]{};
        auto r=solver.solve(epoch,coldStart?zero:reference,systems);
        double enu[3]={NAN,NAN,NAN},vel[3]={NAN,NAN,NAN};
        if (r.positionValid) {
            ++np;
            if (haveReference) {
                double d[3]; for (int i=0; i<3; ++i) d[i]=r.pos[i]-reference[i];
                ecefVectorToEnu(d,reference,enu);
            }
            if (r.velocityValid) { ++nv; ecefVectorToEnu(r.vel,haveReference?reference:r.pos,vel); }
        }
        std::cout<<r.time.week<<','<<std::fixed<<std::setprecision(3)<<r.time.sec<<','<<r.status<<','
                 <<r.nSat<<','<<r.nVel<<','<<(r.positionValid?r.pdop:NAN);
        for (double x:r.pos) std::cout<<','<<(r.positionValid?x:NAN);
        for (double x:enu) std::cout<<','<<x;
        for (double x:vel) std::cout<<','<<x;
        std::cout<<','<<r.nRejected;
        std::cout<<'\n';
    }
    std::cerr<<"Position: "<<np<<'/'<<observations.size()<<"; velocity: "<<nv<<'/'<<observations.size()<<"\n";
    return np;
}
void spp_process(const std::vector<std::vector<ObsData>>& observations,
                 const std::vector<NavData>& navigation, const double* reference) {
    processObservations(observations,navigation,reference,"G");
}
