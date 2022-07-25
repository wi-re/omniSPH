#include <pybind11/pybind11.h>
#include "../../simulation/SPH.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i + j;
}

struct Pet {
    Pet(const std::string &name) : name(name) { }
    void setName(const std::string &name_) { name = name_; }
    const std::string &getName() const { return name; }

    std::string name;
};

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    py::class_<Pet>(m, "Pet")
        .def(py::init<const std::string &>())
        .def("setName", &Pet::setName)
        .def("getName", &Pet::getName);
        
#define fn(x) .def(#x, &SPHSimulation::x)

    py::class_<SPHSimulation>(m,"SPHSimulation")
    .def(py::init<std::string>())
    fn(timestep)
    fn(neighborList)
    fn(fillCells)
    fn(resetFrame)
    fn(density)
    fn(Integrate)
    fn(emitParticles)
    
    fn(computeAlpha)
    fn(computeSourceTerm)
    fn(computeAcceleration)
    fn(updatePressure)
    fn(computeBoundaryTrianglePressure)
    fn(computeBoundaryPressure)
    fn(predictVelocity)
    fn(updateVelocity)
    fn(divergenceSolve)
    fn(densitySolve)

    fn(XSPH)
    fn(computeVorticity)
    fn(refineVorticity)
    fn(externalForces)
    fn(BXSPH)

    fn(dump)

    fn(getScalar)
    fn(getVec)
    fn(getInteger)
    fn(getBoolean)

    fn(setScalar)
    fn(setVec)
    fn(setInteger)
    fn(setBoolean)

    .def_readwrite("fluidPosition", &SPHSimulation::fluidPosition)
    .def_readwrite("fluidVelocity", &SPHSimulation::fluidVelocity)
    .def_readwrite("fluidAccel", &SPHSimulation::fluidAccel)
    .def_readwrite("fluidPredVelocity", &SPHSimulation::fluidPredVelocity)
    .def_readwrite("fluidPredAccel", &SPHSimulation::fluidPredAccel)
    .def_readwrite("fluidPredPosition", &SPHSimulation::fluidPredPosition)

    .def_readwrite("fluidDensity", &SPHSimulation::fluidDensity)
    .def_readwrite("fluidVorticity", &SPHSimulation::fluidVorticity)
    .def_readwrite("fluidAngularVelocity", &SPHSimulation::fluidAngularVelocity)
    .def_readwrite("fluidArea", &SPHSimulation::fluidArea)
    .def_readwrite("fluidRestDensity", &SPHSimulation::fluidRestDensity)
    .def_readwrite("fluidSupport", &SPHSimulation::fluidSupport)

    .def_readwrite("fluidAlpha", &SPHSimulation::fluidAlpha)
    .def_readwrite("fluidActualArea", &SPHSimulation::fluidActualArea)
    .def_readwrite("fluidPressure1", &SPHSimulation::fluidPressure1)
    .def_readwrite("fluidPressure2", &SPHSimulation::fluidPressure2)
    .def_readwrite("fluidBoundaryPressure", &SPHSimulation::fluidBoundaryPressure)
    .def_readwrite("fluidSourceTerm", &SPHSimulation::fluidSourceTerm)
    .def_readwrite("fluidDpDt", &SPHSimulation::fluidDpDt)
    .def_readwrite("fluidSourceTerm", &SPHSimulation::fluidSourceTerm)
    .def_readwrite("fluidDensityStar", &SPHSimulation::fluidDensityStar)
    .def_readwrite("fluidPriorPressure", &SPHSimulation::fluidPriorPressure)

    .def_readwrite("fluidUID", &SPHSimulation::fluidUID)
    
    .def_readonly("fluidNeighbors", &SPHSimulation::fluidNeighborList)
    .def_readonly("boundaryParticles", &SPHSimulation::boundaryParticles);


    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
