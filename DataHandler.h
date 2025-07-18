#pragma once
#include <string>
#include <vector>

struct DataPoint {
    double Q2, xb, t, phi;
    double BSA, BSA_err;
    double XS, XS_err;
};

class DataHandler {
public:
    explicit DataHandler(const std::string& filename);
    const std::vector<DataPoint>& GetData() const;

private:
    std::vector<DataPoint> data;
    void LoadCSV(const std::string& filename);
};
