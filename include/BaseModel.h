#pragma once
#include <string>

class BaseModel {
public:
    virtual ~BaseModel() {}
    virtual void setKinematics(double Q2, double t, double xb, double E) = 0;
    virtual double getBSA() const = 0;
    virtual double getCrossSection() const = 0;
    virtual std::string getName() const = 0;
};
