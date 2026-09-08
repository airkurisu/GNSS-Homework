#include "rinex.h"
#include "spp.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
double hour(const std::string& value) {
    size_t end=0; double h=std::stod(value,&end);
    if (end!=value.size() || !std::isfinite(h) || h<0 || h>24) throw std::runtime_error("Hour must be between 0 and 24");
    return h;
}
}
int main(int argc, char** argv) {
    std::ofstream output;
    auto* original=std::cout.rdbuf();
    try {
        std::string nav="data/brdc1590.24p",obs="data/jfng1590.24o",systems="GCER",out;
        double start=0,end=24;
        bool cold=false;
        for (int i=1; i<argc; ++i) {
            std::string arg=argv[i];
            if (arg=="--help") {
                std::cout<<"Usage: mygnss10 [--nav FILE] [--obs FILE] [--mode gps|multi]\n"
                    <<"                [--systems GCER] [--start HOUR] [--end HOUR]\n"
                    <<"                [--output CSV] [--cold-start]\n"
                    <<"Hours are GPST, start/end inclusive. Default: all epochs, GPS/BDS/Galileo/GLONASS.\n";
                return 0;
            }
            if (arg=="--cold-start") { cold=true; continue; }
            if (i+1>=argc) throw std::runtime_error("Missing value for "+arg);
            std::string value=argv[++i];
            if (arg=="--nav") nav=value;
            else if (arg=="--obs") obs=value;
            else if (arg=="--output") out=value;
            else if (arg=="--start") start=hour(value);
            else if (arg=="--end") end=hour(value);
            else if (arg=="--mode") {
                if (value!="gps" && value!="multi") throw std::runtime_error("Mode must be gps or multi");
                systems=value=="gps"?"G":"GCER";
            } else if (arg=="--systems") {
                if (value.empty() || value.find_first_not_of("GCER")!=std::string::npos)
                    throw std::runtime_error("Systems must contain only G, C, E, R");
                systems=value;
            } else throw std::runtime_error("Unknown option: "+arg);
        }
        if (start>end) throw std::runtime_error("Start hour must not exceed end hour");
        auto navigation=readNavFile(nav);
        double reference[3]{};
        auto observations=readObsFile(obs,reference);
        observations.erase(std::remove_if(observations.begin(),observations.end(),[&](const std::vector<ObsData>& epoch) {
            if (epoch.empty()) return true;
            double h=std::fmod(epoch[0].time.sec,86400.0)/3600.0;
            return h<start-1e-10 || h>end+1e-10;
        }),observations.end());
        if (navigation.empty() || observations.empty()) throw std::runtime_error("No navigation/observation data in selected interval");
        if (!out.empty()) {
            output.open(out);
            if (!output) throw std::runtime_error("Cannot open output: "+out);
            std::cout.rdbuf(output.rdbuf());
        }
        const int solved=processObservations(observations,navigation,reference,systems,cold);
        std::cout.flush();
        bool ok=bool(std::cout);
        std::cout.rdbuf(original);
        if (!ok) throw std::runtime_error("Failed to write output");
        return solved>0?0:2;
    } catch (const std::exception& e) {
        std::cout.rdbuf(original);
        std::cerr<<"Error: "<<e.what()<<'\n';
        return 1;
    }
}
