#include "DataHandler.h"
#include "Fitter.h"

int main() {
    DataHandler dh("data/sample_data.csv");
    Fitter fitter(dh.GetData());
    fitter.Fit();
    return 0;
}
