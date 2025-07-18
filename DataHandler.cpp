#include "DataHandler.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

DataHandler::DataHandler(const std::string& filename) {
    LoadCSV(filename);
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
