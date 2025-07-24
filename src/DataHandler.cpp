#include "DataHandler.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

DataHandler::DataHandler() {
}

const std::vector<DataPoint>& DataHandler::GetData() const {
    return data;
}

void DataHandler::LoadCSV(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        DataPoint dp;
        char comma;
        ss >> dp.Q2 >> comma >> dp.xb >> comma >> dp.t >> comma >> dp.phi >> comma
           >> dp.BSA >> comma >> dp.BSA_err >> comma >> dp.XS >> comma >> dp.XS_err;
        data.push_back(dp);
    }
}
std::vector<KinematicBlock> DataHandler::loadData(const std::string& filename) {
    std::ifstream infile(filename);
    std::string line;

    std::map<std::tuple<double, double, double>, std::vector<InputDataPoint>> grouped;
     std::cout<<" reading file" <<std::endl;
    while (std::getline(infile, line)) {
         //std::cout<<" reading lines" <<std::endl;
        if (line.empty() || line[0] == '#') {
           std::cout<<" line[0]" << line[0]<<std::endl;
           continue;
        }

        std::istringstream iss(line);
        double xB, Q2, t, phi;
        double xlu, xlu_stat, xlu_syst;
        double cs, cs_stat, cs_syst;

        if (!(iss >> xB >> Q2 >> t >> phi >> xlu >> xlu_stat >> xlu_syst >> cs >> cs_stat >> cs_syst)){
            continue; // skip if not enough values
        }
        // Optional: skip zero-error placeholders
        if ((xlu_stat == 0 && xlu_syst == 0) && (cs_stat == 0 && cs_syst == 0)){
            std::cout<<" phi-3" <<phi<<std::endl;
            continue;
        }

        InputDataPoint point;
        point.phi     = phi;
        point.bsa     = xlu;
        point.bsa_err = std::sqrt(xlu_stat * xlu_stat + xlu_syst * xlu_syst);
        point.cs      = cs;
        point.cs_err  = std::sqrt(cs_stat * cs_stat + cs_syst * cs_syst);

        grouped[{xB, Q2, -t}].push_back(point);  // Negate t to convert from -t to t

        //std::cout<<" "<<xB <<" "<< Q2 <<" "<< t <<" "<<phi <<" "<< xlu <<" "<< xlu_stat <<" "<< xlu_syst <<" "<< cs <<" "<< cs_stat <<" "<< cs_syst<<std::endl;
    }

    std::vector<KinematicBlock> blocks;
    for (const auto& [key, points] : grouped) {
        const auto& [xB, Q2, t] = key;
        blocks.push_back({xB, Q2, t, points});
    }

    return blocks;
}

