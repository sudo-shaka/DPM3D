/*
 * =====================================================================================
 *
 *       Filename:  ECM_Wrapper.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  09/01/2022 02:08:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */
#include "../include/ECM.hpp"
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_ECM(py::module &m){
  py::class_<DPM3D::ECM>(m, "ECM")
    .def(py::init<std::vector<glm::dvec3>>(),
         py::arg("Positions")
         )
    .def(py::init<double,double,int>(),
         py::arg("Z"),
         py::arg("L"),
         py::arg("NP")
         )
    .def_readwrite("Positions",&DPM3D::ECM::Positions)
    .def_readonly("Forces",&DPM3D::ECM::Forces)
    .def_readonly("NP",&DPM3D::ECM::NP)
    ;
}
