#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Cell(py::module &);
void init_Tissue(py::module &);
void init_dvec3(py::module &);

namespace dpm {

    PYBIND11_MODULE(DPM3D, m){
        m.doc() = "Deformable Particle Model";
        init_Cell(m);
        init_Tissue(m);
        init_dvec3(m);
    }
}
