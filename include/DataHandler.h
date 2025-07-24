#pragma once
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

struct InputDataPoint {
    double phi;
    double bsa;
    double cs;
    double bsa_err;
    double cs_err;
};

struct DataPoint {
    double Q2, xb, t, phi;
    double BSA, BSA_err;
    double XS, XS_err;
};

struct KinematicBlock {
    double xB, Q2, t;
    std::vector<InputDataPoint> points;
};
class DataHandler {
public:
    explicit DataHandler();
    const std::vector<DataPoint>& GetData() const;
    std::vector<KinematicBlock> loadData(const std::string& filename);

private:
    std::vector<DataPoint> data;
    void LoadCSV(const std::string& filename);
};
