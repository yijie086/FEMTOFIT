#pragma once
#include "BaseModel.h"

class KM15Model : public BaseModel {
public:
    KM15Model();
    ~KM15Model();

    void setKinematics(double Q2, double t, double xb, double E) override;
    double getBSA() const override;
    double getCrossSection() const override;
    std::string getName() const override;

private:
    double lastBSA = 0;
    double lastXS = 0;
};
