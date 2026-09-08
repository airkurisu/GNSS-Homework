#include "rinex.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
std::string field(const std::string& s, size_t start, size_t n) {
    return start<s.size() ? s.substr(start,n) : "";
}
bool blank(const std::string& s) { return s.find_first_not_of(" \r\t")==std::string::npos; }
double number(std::string s, double missing=0) {
    if (blank(s)) return missing;
    for (char& c:s) if (c=='D' || c=='d') c='E';
    size_t end=0;
    double value=std::stod(s,&end);
    if (!blank(s.substr(end)) || !std::isfinite(value)) throw std::runtime_error("Invalid numeric field: "+s);
    return value;
}
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void version(std::ifstream& file, char type) {
    std::string line;
    require(bool(std::getline(file,line)),"Empty RINEX file");
    double v=number(field(line,0,9));
    require(v>=3.0 && v<=3.04+1e-6,"Only RINEX 3.00-3.04 is supported");
    require(field(line,20,1)==std::string(1,type),"Wrong RINEX file type");
}
int satellite(const std::string& line) {
    require(line.size()>=3 && line[1]>='0' && line[1]<='9' && line[2]>='0' && line[2]<='9',"Invalid satellite record");
    int prn=(line[1]-'0')*10+line[2]-'0';
    require(prn>0,"Invalid PRN");
    return prn;
}
// Deliberately restrict to signal pairs whose broadcast clock convention is modeled.
const std::map<char,std::vector<std::string>> firstCodes={
    {'G',{"C1W","C1C","C1P","C1X","C1S","C1L"}},
    {'E',{"C1X","C1C","C1B"}}, {'C',{"C2I","C2Q","C2X"}},
    {'R',{"C1P","C1C"}}
};
const std::map<char,std::vector<std::string>> secondCodes={
    {'G',{"C2W","C2P","C2X","C2L","C2S","C2C"}},
    {'E',{"C5X","C5Q","C5I","C7X","C7Q","C7I"}},
    {'C',{"C6I","C6Q","C6X","C7I","C7Q","C7X"}}, {'R',{"C2P","C2C"}}
};
void readTypes(std::ifstream& file, std::string line, std::map<char,std::vector<std::string>>& types) {
    const char sys=line[0];
    int count=static_cast<int>(number(field(line,3,3)));
    require(count>0 && count<=128,"Invalid observation type count");
    std::vector<std::string> codes;
    for (int i=0; i<count; ++i) {
        if (i>0 && i%13==0) {
            require(bool(std::getline(file,line)) && line.find("SYS / # / OBS TYPES")!=std::string::npos,
                    "Missing observation type continuation");
        }
        std::string code=field(line,7+4*(i%13),3);
        require(code.size()==3 && !blank(code),"Missing observation type");
        codes.push_back(code);
    }
    types[sys]=codes;
}
}

std::vector<NavData> readNavFile(const std::string& filename) {
    std::ifstream file(filename);
    require(file.is_open(),"Cannot open navigation file: "+filename);
    version(file,'N');
    std::string line;
    bool ended=false;
    while (std::getline(file,line)) if (line.find("END OF HEADER")!=std::string::npos) { ended=true; break; }
    require(ended,"Missing navigation END OF HEADER");
    std::vector<NavData> result;
    while (std::getline(file,line)) {
        if (blank(line)) continue;
        NavData nav;
        nav.sat=satellite(line); nav.sys=line[0];
        require(std::string("GRECJSI").find(nav.sys)!=std::string::npos,"Unknown navigation system");
        double ep[6]{};
        std::istringstream epoch(field(line,4,19));
        for (double& x:ep) require(bool(epoch>>x),"Invalid navigation epoch");
        nav.toc=epoch2time(ep);
        nav.a0=number(field(line,23,19),NAN);
        nav.a1=number(field(line,42,19),NAN);
        nav.a2=number(field(line,61,19),NAN);
        require(std::isfinite(nav.a0+nav.a1+nav.a2),"Incomplete navigation clock record");
        const int rows=(nav.sys=='R' || nav.sys=='S')?3:7;
        double b[7][4]{};
        for (int i=0; i<rows; ++i) {
            require(bool(std::getline(file,line)) && line.size()>=4 && blank(line.substr(0,4)),"Truncated navigation record");
            for (int j=0; j<4; ++j) b[i][j]=number(field(line,4+19*j,19));
            // The final row contains optional/spare fields; earlier rows must retain core fields.
            if (i<rows-1) require(line.size()>=61,"Truncated navigation data row");
        }
        if (nav.sys=='R') {
            nav.toc=utc2gpst(nav.toc); nav.toe=nav.toc;
            for (int j=0; j<3; ++j) {
                nav.gloPos[j]=b[j][0]*1000; nav.gloVel[j]=b[j][1]*1000; nav.gloAcc[j]=b[j][2]*1000;
            }
            nav.health=static_cast<int>(b[0][3]); nav.freq_num=static_cast<int>(b[1][3]);
            require(nav.freq_num>=-7 && nav.freq_num<=6,"Invalid GLONASS frequency channel");
        } else if (nav.sys=='G' || nav.sys=='E' || nav.sys=='C') {
            nav.crs=b[0][1]; nav.delta_n=b[0][2]; nav.M0=b[0][3];
            nav.cuc=b[1][0]; nav.e=b[1][1]; nav.cus=b[1][2]; nav.sqrtA=b[1][3];
            nav.toeSow=b[2][0]; nav.cic=b[2][1]; nav.Omega0=b[2][2]; nav.cis=b[2][3];
            nav.i0=b[3][0]; nav.crc=b[3][1]; nav.omega=b[3][2]; nav.OmegaDot=b[3][3];
            nav.IDot=b[4][0]; nav.dataSources=static_cast<int>(b[4][1]);
            nav.health=static_cast<int>(b[5][1]); nav.tgd=b[5][2]; nav.tgd2=b[5][3];
            nav.toe={static_cast<int>(b[4][2]),nav.toeSow};
            if (nav.sys=='C') {
                nav.toe.week+=1356; // BDT epoch 2006-01-01
                nav.toe=timeadd(nav.toe,14); nav.toc=timeadd(nav.toc,14);
            }
            require(nav.sqrtA>0 && nav.e>=0 && nav.e<1 && nav.toeSow>=0 && nav.toeSow<604800,
                    "Invalid Kepler ephemeris");
        } else continue; // consume complete unsupported records before skipping them
        result.push_back(nav);
    }
    std::cerr<<"Loaded "<<result.size()<<" navigation records.\n";
    return result;
}

std::vector<std::vector<ObsData>> readObsFile(const std::string& filename, double* approxPos) {
    std::ifstream file(filename);
    require(file.is_open(),"Cannot open observation file: "+filename);
    version(file,'O');
    std::fill(approxPos,approxPos+3,0.0);
    std::map<char,std::vector<std::string>> types;
    std::string line, timeSystem="GPS";
    bool ended=false;
    while (std::getline(file,line)) {
        if (line.find("APPROX POSITION XYZ")!=std::string::npos)
            for (int i=0; i<3; ++i) approxPos[i]=number(field(line,14*i,14));
        if (line.find("SYS / # / OBS TYPES")!=std::string::npos) readTypes(file,line,types);
        if (line.find("TIME OF FIRST OBS")!=std::string::npos) {
            std::string sys=field(line,48,3);
            if (!blank(sys)) timeSystem=sys;
        }
        if (line.find("SYS / SCALE FACTOR")!=std::string::npos)
            require(number(field(line,2,4))==1,"Non-unit observation scale factors are unsupported");
        if (line.find("END OF HEADER")!=std::string::npos) { ended=true; break; }
    }
    require(ended && !types.empty(),"Missing observation header/types");
    require(timeSystem=="GPS" || timeSystem=="GAL" || timeSystem=="BDT" || timeSystem=="GLO" || timeSystem=="UTC",
            "Unsupported observation time system: "+timeSystem);
    std::vector<std::vector<ObsData>> result;
    GPSTime previous{}; bool havePrevious=false;
    while (std::getline(file,line)) {
        if (blank(line)) continue;
        require(line[0]=='>',"Expected observation epoch");
        std::istringstream header(line.substr(1));
        double ep[6]{};
        for (double& x:ep) require(bool(header>>x),"Invalid observation epoch");
        int flag=0, count=0;
        require(bool(header>>flag>>count) && count>=0 && count<=999,"Invalid epoch flag/count");
        require(flag>=0 && flag<=6,"Unknown epoch event flag");
        require(flag!=4,"Mid-file header changes are unsupported");
        if (flag>=2 && flag<=5) {
            for (int i=0; i<count; ++i) require(bool(std::getline(file,line)),"Truncated event record");
            continue;
        }
        double clockOffset=0;
        if (header>>clockOffset) require(clockOffset==0,"Nonzero epoch receiver clock offsets are unsupported");
        GPSTime time=epoch2time(ep);
        if (timeSystem=="BDT") time=timeadd(time,14);
        if (timeSystem=="GLO" || timeSystem=="UTC") time=utc2gpst(time);
        if (flag!=6) {
            require(!havePrevious || timediff(time,previous)>0,"Observation epochs must increase");
            previous=time; havePrevious=true;
        }
        std::vector<ObsData> epoch;
        std::set<std::string> seen;
        for (int i=0; i<count; ++i) {
            require(bool(std::getline(file,line)),"Truncated observation epoch");
            ObsData obs; obs.time=time; obs.sat=satellite(line); obs.sys=line[0];
            require(seen.insert(line.substr(0,3)).second,"Duplicate satellite in epoch");
            auto type=types.find(obs.sys);
            require(type!=types.end(),"Missing observation types for satellite");
            const auto& codes=type->second;
            if (flag==6 || firstCodes.count(obs.sys)==0) continue;
            // RINEX 3 stores one satellite per (possibly very long) line; omitted trailing fields are blank.
            auto value=[&](const std::string& code, double missing) {
                auto it=std::find(codes.begin(),codes.end(),code);
                return it==codes.end()?missing:number(field(line,3+16*static_cast<size_t>(it-codes.begin()),14),missing);
            };
            auto choose=[&](const std::vector<std::string>& priorities, std::string& code, double& p, double& doppler) {
                for (const auto& candidate:priorities) {
                    double v=value(candidate,0);
                    if (v>0) { code=candidate; p=v; doppler=value("D"+candidate.substr(1),NAN); return; }
                }
            };
            choose(firstCodes.at(obs.sys),obs.codeP1,obs.P1,obs.D1);
            choose(secondCodes.at(obs.sys),obs.codeP2,obs.P2,obs.D2);
            if (obs.P1>0 && obs.P2>0) epoch.push_back(obs);
        }
        if (flag!=6) {
            // A time-only sentinel preserves the epoch for explicit failure reporting.
            if (epoch.empty()) { ObsData missing; missing.time=time; epoch.push_back(missing); }
            result.push_back(epoch);
        }
    }
    std::cerr<<"Loaded "<<result.size()<<" epochs.\n";
    return result;
}
