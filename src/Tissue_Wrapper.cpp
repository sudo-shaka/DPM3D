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
    .def_readwrite("L",&DPM3D::Tissue::L)
    .def_readonly("phi0",&DPM3D::Tissue::phi0)
    .def_readwrite("Kc",&DPM3D::Tissue::Kc)
    .def_readonly("NCELLS",&DPM3D::Tissue::NCELLS)
    .def_readonly("Cells",&DPM3D::Tissue::Cells)
    .def_readonly("overlaps",&DPM3D::Tissue::overlaps)
    .def("InteractingForceUpdate",&DPM3D::Tissue::RetractingForceUpdate)
    .def("UpdateShapeForces",&DPM3D::Tissue::UpdateShapeForces)
    .def("EulerUpdate",py::overload_cast<double>(&DPM3D::Tissue::EulerUpdate),py::arg("dt"))
    .def("EulerUpdate",py::overload_cast<int,double>(&DPM3D::Tissue::EulerUpdate),py::arg("nsteps"),py::arg("dt"))
    .def("FindOverlaps",&DPM3D::Tissue::FindOverlaps)
    .def("MonolayerDisperse",&DPM3D::Tissue::MonolayerDisperse)
    .def("TissueDisperse",&DPM3D::Tissue::TissueDisperse)
    ;
}
