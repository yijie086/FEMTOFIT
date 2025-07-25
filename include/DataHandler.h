/*#pragma once
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

struct FitResult {
    double xB;
    double t;
    double xi;  // Optional, calculated from xB and Q2
    double ImH;
    double ReH;
    double ImHerr;
    double ReHerr;
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
*/
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
    double Q2;
    double xb;
    double t;
    double phi;
    double BSA;
    double BSA_err;
    double XS;
    double XS_err;
};

struct KinematicBlock {
    double xB;
    double Q2;
    double t;
    std::vector<InputDataPoint> points;
};

struct FitResult {
    double xB;
    double t;
    double xi;      // Computed from xB, Q2
    double ImH;
    double ReH;
    double ImHerr;
    double ReHerr;
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