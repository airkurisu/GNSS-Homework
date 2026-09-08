// Read "week tow satellite" queries on stdin; emit broadcast states as CSV.
// Evaluate at the query epoch itself (not a receive/transmit-time transformation).
#include "spp.h"
#include "orbit.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc!=2) throw std::runtime_error("Usage: gnss_state_dump NAV < queries.txt > states.csv");
        SPPSolver solver(readNavFile(argv[1]));
        std::cout<<"sat,week,tow,status,x,y,z,vx,vy,vz,clock,drift,clock_polynomial,drift_polynomial,toe_week,toe_sow,toc_week,toc_sow,tgd,tgd2,data_sources\n";
        int week; double tow; std::string sat;
        while (std::cin>>week>>tow>>sat) {
            if (sat.size()!=3) throw std::runtime_error("Invalid query satellite");
            ObsData o; o.sys=sat[0]; o.sat=std::stoi(sat.substr(1)); o.time={week,tow};
            o.codeP1=o.sys=='C'?"C2I":"C1P";
            o.codeP2=o.sys=='C'?"C6I":(o.sys=='E'?"C5Q":"C2P");
            const NavData* nav=solver.selectEphemeris(o);
            std::cout<<sat<<','<<week<<','<<std::setprecision(15)<<tow;
            if (!nav) { std::cout<<",no_ephemeris,,,,,,,,,,,,,,,,,\n"; continue; }
            double p[3],v[3],clock,drift;
            satPosVel(o.time,*nav,p,v,&clock,&drift);
            double tc=timediff(o.time,nav->toc);
            double poly=nav->a0+nav->a1*tc+(o.sys=='R'?0:nav->a2*tc*tc);
            double polyDrift=nav->a1+(o.sys=='R'?0:2*nav->a2*tc);
            std::cout<<",ok";
            for (double x:p) std::cout<<','<<x;
            for (double x:v) std::cout<<','<<x;
            std::cout<<','<<clock<<','<<drift<<','<<poly<<','<<polyDrift<<','<<nav->toe.week<<','
                     <<nav->toe.sec<<','<<nav->toc.week<<','<<nav->toc.sec<<','<<nav->tgd<<','
                     <<nav->tgd2<<','<<nav->dataSources<<'\n';
        }
        return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
