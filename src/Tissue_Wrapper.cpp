#include "../include/Tissue.hpp"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py=pybind11;

void init_Tissue(py::module &m){
    py::class_<DPM3D::Tissue>(m,"Tissue")
    .def(py::init<std::vector<DPM3D::Cell>,double>(),
    py::arg("Cells"),
    py::arg("phi0")
    )
    .def_readonly("phi0",&DPM3D::Tissue::phi0)
    .def_readwrite("Kc",&DPM3D::Tissue::Kc)
    .def_readonly("NCELLS",&DPM3D::Tissue::NCELLS)
    .def_readonly("Cells",&DPM3D::Tissue::Cells)
    .def("InteractingForceUpdate",&DPM3D::Tissue::InteractingForceUpdate)
    ;
}