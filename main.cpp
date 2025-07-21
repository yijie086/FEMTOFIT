//STL
#include <iostream>

#include "DataHandler.h"
#include "Fitter.h"

/// models for fitting
#include "Models.h"
#include "KM15Model.h"

int main() {
    //DataHandler dh("data/sample_data.csv");
    //Fitter fitter(dh.GetData());
    //fitter.Fit();

    Models modelManager;
    modelManager.registerModel(std::make_shared<KM15Model>());

    auto model = modelManager.getModel("KM15Model");
    model->setKinematics(2.0, -0.3, 0.1, 10.6);

    std::cout << "BSA: " << model->getBSA() << std::endl;
    std::cout << "Cross Section: " << model->getCrossSection() << std::endl;

    return 0;
}
