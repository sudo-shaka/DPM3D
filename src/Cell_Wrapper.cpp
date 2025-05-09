#include "Cell3D.hpp"
#include "ECM.hpp"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_dvec3(py::module &m){
    py::class_<glm::dvec3>(m,"dvec3")
    .def(py::init<>())
    .def(py::init<double,double,double>(),py::arg("x"),py::arg("y"),py::arg("z"))
    .def_readwrite("x",&glm::dvec3::x)
    .def_readwrite("y",&glm::dvec3::y)
    .def_readwrite("z",&glm::dvec3::z)
    ;
}

void init_Cell(py::module &m){
    py::class_<DPM3D::Cell>(m, "Cell")
    .def(py::init<double, double, double,double, int,double, double, double,double>(),
    py::arg("x1"),
    py::arg("y1"),
    py::arg("z1"),
    py::arg("calA0"),
    py::arg("VertexRecursion"),
    py::arg("r0"),
    py::arg("Kv"),
    py::arg("Ka"),
    py::arg("Kb")
    )
    .def(py::init<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<std::vector<int>>>(),
      py::arg("X"),
      py::arg("Y"),
      py::arg("Z"),
      py::arg("FaceIdx")
    )
    .def_readonly("NV",&DPM3D::Cell::NV)
    .def_readonly("isJunction",&DPM3D::Cell::isJunction)
    .def_readonly("isFocalAdhesion",&DPM3D::Cell::isFocalAdh)
    .def_readwrite("a0",&DPM3D::Cell::a0)
    .def_readwrite("v0",&DPM3D::Cell::v0)
    .def_readwrite("Kb",&DPM3D::Cell::Kb)
    .def_readwrite("Ka",&DPM3D::Cell::Ka)
    .def_readwrite("Kl",&DPM3D::Cell::Kl)
    .def_readwrite("Kv",&DPM3D::Cell::Kb)
    .def_readwrite("Positions",&DPM3D::Cell::Positions)
    .def_readwrite("TriangleIndex",&DPM3D::Cell::FaceIndices)
    .def_readwrite("Velocities",&DPM3D::Cell::Velocities)
    .def_readwrite("Forces",&DPM3D::Cell::Forces)
    .def_readwrite("Ks",&DPM3D::Cell::Ks)
    .def("ResetForces",&DPM3D::Cell::ResetForces)
    .def("UpdateShapeForces",&DPM3D::Cell::ShapeForceUpdate)
    .def("UpdateSurfaceAreaForce",&DPM3D::Cell::AreaForceUpdate)
    .def("UpdateVolumeForce",&DPM3D::Cell::VolumeForceUpdate)
    .def("UpdateBendingForce",&DPM3D::Cell::BendingForceUpdate)
    .def("GetShapeForces",&DPM3D::Cell::GetShapeForces)
    .def("GetVolumeForces",&DPM3D::Cell::GetVolumeForces)
    .def("GetAreaForces",&DPM3D::Cell::GetAreaForces)
    .def("GetPositions",&DPM3D::Cell::GetPositions)
    .def("GetForces",&DPM3D::Cell::GetForces)
    .def("GetBendingForces",&DPM3D::Cell::GetBendingForces)
    .def("FindJunctions",&DPM3D::Cell::FindJunctions)
    .def("GetJunctionPositions",&DPM3D::Cell::getJunctionPoints)
    .def("VertexFIRE",py::overload_cast<double, double, int, double>(&DPM3D::Cell::FIREMinimization),py::arg("alpha0"),py::arg("dt"),py::arg("itmax"),py::arg("Ftol"))
    .def("EulerUpdate",py::overload_cast<double>(&DPM3D::Cell::EulerUpdate),py::arg("dt"))
    .def("EulerUpdate",py::overload_cast<int,double>(&DPM3D::Cell::EulerUpdate),py::arg("nsteps"),py::arg("dt"))
    .def("StickToSurface",py::overload_cast<double,double>(&DPM3D::Cell::StickToSurface),py::arg("minz"),py::arg("mindist"))
    .def("StickToSurface",py::overload_cast<DPM3D::ECM,double>(&DPM3D::Cell::StickToSurface),py::arg("ECM"),py::arg("mindist"))
    .def("StickToSurfaceCatch",py::overload_cast<double,double>(&DPM3D::Cell::StickToSurfaceCatch),py::arg("minz"),py::arg("mindistance"))
    .def("StickToSurfaceCatch",py::overload_cast<DPM3D::ECM,double>(&DPM3D::Cell::StickToSurfaceCatch),py::arg("ECM"),py::arg("mindist"))
    .def("StickToSurfaceSlip",py::overload_cast<double,double>(&DPM3D::Cell::StickToSurfaceSlip),py::arg("minz"),py::arg("mindist"))
    .def("StickToSurfaceSlip",py::overload_cast<DPM3D::ECM,double>(&DPM3D::Cell::StickToSurfaceSlip),py::arg("ECM"),py::arg("mindist"))
    .def("SurfaceGradient",&DPM3D::Cell::SurfaceGradient)
    .def("RepelSurface",&DPM3D::Cell::RepelSurface)
    .def("ExtendVertex",py::overload_cast<int,double>(&DPM3D::Cell::ExtendVertex),py::arg("vertexidx"),py::arg("Force"))
    .def("GetVolume",&DPM3D::Cell::GetVolume)
    .def("GetSA",&DPM3D::Cell::GetSurfaceArea)
    .def("GetCalA0",&DPM3D::Cell::GetCalA0)
    .def("PointInside",&DPM3D::Cell::pointInside)
    ;
}
