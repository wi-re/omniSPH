#pragma once
//#define _USE_MATH_DEFINES
#include <array>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <tools/ParameterManager.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>

#include <yaml-cpp/yaml.h>

// comment this line out to use double precision for everything in the simulation
//#define USE_SINGLE
#ifndef USE_SINGLE
using vec = Eigen::Vector2d;
using scalar = double;
using matrix = Eigen::Matrix2d;
using complex = std::complex<scalar>;
#else
using vec = Eigen::Vector2f;
using scalar = float;
using matrix = Eigen::Matrix2f;
using complex = std::complex<scalar>;
#endif
// using define for easier usage of time
using clk = std::chrono::high_resolution_clock;
// constexpr inline scalar domainScale = 100.0;

constexpr inline double double_pi = 3.141592653589793238462643383279502884L;

inline std::ofstream summaryFile;
inline bool summaryFileOpen = false;
//// kernel Parameters
constexpr inline int32_t targetNeighbors = 20;
constexpr inline scalar kernelSize = (scalar)1.778002;
//// numerical Parameters
constexpr inline scalar epsilon((scalar)1e-7);

struct Triangle {
    vec v0 = vec(0, 0), v1 = vec(0, 0), v2 = vec(0, 0);
};
std::tuple<bool, vec, scalar, scalar, vec> interactTriangle(vec p, scalar support, Triangle tri);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentric(vec p, scalar support, scalar rho0, Triangle tri, scalar rhoi, scalar fi, scalar f0, scalar f1, scalar f2);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentricSimple(vec p, scalar support, scalar rho0, Triangle tri,scalar f0, scalar f1, scalar f2);
enum struct shape_t{
    rectangular, spherical
};
enum struct emitter_t{
    oneTime, inlet, outlet, velocitySource
};

struct fluidSource{
    vec emitterMin;
    vec emitterMax;

    scalar emitterDensity = 998.0;
    scalar emitterRadius;

    scalar emitterOffset = 0.;
    scalar timeLimit = -1.;
    scalar emitterRampTime = -1.;

    vec emitterVelocity = vec(0,0);

    shape_t shape = shape_t::rectangular;
    emitter_t emitter = emitter_t::oneTime;

    std::vector<vec> genParticles() const;
};

struct gravitySource{
    bool pointSource = false;
    vec direction = vec(0,-1.);
    vec location = vec(0.,0.);
    scalar magnitude = -9.8;
};

class SPHSimulation{
    fs::path actualPath;
    bool initDump = false;
    int32_t exportFrameCounter = 0;
    int32_t lastExport = 0.;

    std::size_t cellsX, cellsY;
    vec domainMin, domainMax;
    scalar baseSupport, baseRadius;


    YAML::Node config;
    void initializeSPH();
    void initializeParameters();

    std::vector<vec> genParticles(vec minCoord, vec maxCoord, scalar radius, scalar packing = 0.39960767069134208);
    std::vector<vec> genParticlesCircle(vec minCoord, vec maxCoord, scalar radius, scalar packing = 0.39960767069134208);
    template<typename Func> auto boundaryFunc(const int32_t ptclIndex, Func&& c) {
    auto p = fluidPosition[ptclIndex];
    auto h = fluidSupport[ptclIndex];
    for (auto ti : fluidTriangleNeighborList[ptclIndex]) {
        const auto& tri = boundaryTriangles[ti];
        auto [hit, pb, d, k, gk] = interactTriangle(p, h, tri);
        if (hit)
            c(pb, d, k, gk, true);
    }
}


    std::vector<int32_t>& getCell(scalar x, scalar y);
    std::vector<int32_t>& getBoundaryCell(int32_t xi, int32_t yi);
    std::vector<int32_t>& getTriangleCell(int32_t xi, int32_t yi);

public:
    SPHSimulation(std::string _config = ""){
        if(_config != ""){
            config = YAML::Load(_config);
            initializeParameters();
            initializeSPH();
        }
    }

    std::vector<vec> fluidPosition, fluidVelocity, fluidAccel, fluidPredVelocity, fluidPredAccel, fluidPredPosition;
    std::vector<scalar> fluidDensity, fluidVorticity, fluidAngularVelocity, fluidArea, fluidRestDensity, fluidSupport;
    std::vector<scalar> fluidAlpha, fluidActualArea, fluidPressure1, fluidPressure2, fluidBoundaryPressure, fluidSourceTerm, fluidDpDt, fluidDensityStar, fluidPriorPressure;
    std::vector<std::vector<int32_t>> fluidNeighborList, fluidTriangleNeighborList;
    std::vector<int32_t> fluidUID;
    int32_t fluidCounter = 0;
    std::vector<std::tuple<scalar, scalar, scalar>> boundaryBarycentricPressure;
    std::vector<vec> boundaryPolygon;
    std::vector<std::vector<vec>> boundaryObstacles;

    std::vector<std::vector<int32_t>> cellArray, cellBoundaryArray, cellTriangleArray;

    std::vector<std::pair<vec,vec>> boundaryParticles;

    std::vector<Triangle> boundaryTriangles;

    std::vector<fluidSource> fluidSources;
    std::vector<gravitySource> gravitySources;


    enum struct property_t{
        position, velocity, accel, predVelocity, predAccel, predPosition,
        density, vorticity, angularVelocity, area, restDensity, support,
        alpha, actualArea, pressure1, pressure2, boundaryPressure, sourceTerm, dpdt, rhoStar, priorPressure,
        UID, neighbors
    };

    std::pair<int32_t, int32_t> getCellIdx(scalar x, scalar y);
    std::vector<int32_t>& getCell(int32_t xi, int32_t yi);

    inline std::vector<vec>& getPositions(){return fluidPosition;}

    void timestep();

    void neighborList();
    void fillCells();
    void resetFrame();
    void density();
    void Integrate();
    void emitParticles();

    void computeAlpha(bool density = true);
    void computeSourceTerm(bool density = true);
    void computeAcceleration(bool density = true);
    void updatePressure(bool density = true);
    void computeBoundaryTrianglePressure(bool density = true);
    scalar calculateBoundaryPressureMLS(int32_t i, vec pb, bool density = true);
    void computeBoundaryPressure(bool density = true);
    void predictVelocity(bool density = true);
    void updateVelocity(bool density = true);
    int32_t divergenceSolve();
    int32_t densitySolve();

    void XSPH();
    void computeVorticity();
    void refineVorticity();
    void externalForces();
    void BXSPH();

    void dump();

    std::tuple<scalar, scalar, std::vector<scalar>> colorMap(property_t prop, bool autoMinMax = true, scalar min = 0., scalar max = 1.);
   
    ParameterManager pm;

    fs::path resolveFile(std::string fileName, std::vector<std::string> search_paths = {});
};
inline SPHSimulation simulationState;
