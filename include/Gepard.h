#pragma once
#include <pybind11/embed.h>
namespace py = pybind11;

class Gepard {
public:
    static void call_KM15(double Q2, double t, double xb, double E, double& bsa, double& xs) {
        static Gepard instance;
        py::object pt = instance.gepard.attr("DataPoint")(
            py::arg("xB") = xb,
            py::arg("t") = t,
            py::arg("Q2") = Q2,
            py::arg("phi") = 0,
            py::arg("process") = "ep2epgamma",
            py::arg("exptype") = "fixed target",
            py::arg("in1energy") = E,
            py::arg("in1charge") = -1,
            py::arg("in1polarization") = -1,
            py::arg("observable") = "XS"
        );

        pt.attr("prepare")();
        bsa = instance.km15.attr("ALUI")(pt).cast<double>();
        xs  = instance.km15.attr("XS")(pt).cast<double>();
    }

private:
    Gepard()
        : guard{}
    {
        // Append Gepard's src directory to PYTHONPATH
        py::module_::import("sys").attr("path").attr("insert")(0,
            "/w/hallb-scshelf2102/clas12/singh/Softwares/Gepard/gepard/src");

        gepard = py::module_::import("gepard");
        km15 = py::module_::import("gepard.fits.th_KM15");
    }

    py::scoped_interpreter guard;
    py::module_ gepard;
    py::module_ km15;
};
