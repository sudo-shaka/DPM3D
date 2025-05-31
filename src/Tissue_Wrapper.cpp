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
    .def_readwrite("PBC",&DPM3D::Tissue::PBC)
    .def_readonly("phi0",&DPM3D::Tissue::phi0)
    .def_readwrite("Kre",&DPM3D::Tissue::Kre)
    .def_readwrite("Kat",&DPM3D::Tissue::Kat)
    .def_readwrite("f0",&DPM3D::Tissue::f0)
    .def_readwrite("AttractionMethod",&DPM3D::Tissue::attactionMethod)
    .def_readonly("NCELLS",&DPM3D::Tissue::NCELLS)
    .def_readonly("Cells",&DPM3D::Tissue::Cells)
    .def("UpdateShapeForces",&DPM3D::Tissue::UpdateShapeForces)
    .def("EulerUpdate",py::overload_cast<double>(&DPM3D::Tissue::EulerUpdate),py::arg("dt"))
    .def("EulerUpdate",py::overload_cast<int,double>(&DPM3D::Tissue::EulerUpdate),py::arg("nsteps"),py::arg("dt"))
    .def("FindOverlaps",&DPM3D::Tissue::FindOverlaps)
    .def("MonolayerDisperse",&DPM3D::Tissue::MonolayerDisperse)
    .def("TissueDisperse",&DPM3D::Tissue::TissueDisperse)
    .def("StickToSurface",&DPM3D::Tissue::StickToSurface)
    .def("StickToSurfaceStrech",&DPM3D::Tissue::StickToSurfaceStretch)
    .def("InteractingUpdate",&DPM3D::Tissue::InteractingUpdate)
    .def("InteractECM",&DPM3D::Tissue::InteractECM)
    .def("JunctionSlipForceUpdate",&DPM3D::Tissue::JunctionSlipForceUpdate)
    .def("JunctionCatchForceUpdate",&DPM3D::Tissue::JunctionSlipForceUpdate)
    .def("UpdateJunctions",&DPM3D::Tissue::UpdateJunctions)
    .def("GetVesselPositions",&DPM3D::Tissue::GetVesselPosition)
    .def("applyShearStress",&DPM3D::Tissue::applyShearStress)
    .def("ExportForces",&DPM3D::Tissue::ExportForces)
    .def("ExportPositions",&DPM3D::Tissue::ExportPositions)
    ;
}
