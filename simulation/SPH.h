#pragma once
//#define _USE_MATH_DEFINES
#include <array>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <boost/atomic.hpp>
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
// constexpr inline scalar packing_2D = 0.39960767069134208;
constexpr inline scalar packing_2D = 0.399200743165053487;

struct Triangle {
    vec v0 = vec(0, 0), v1 = vec(0, 0), v2 = vec(0, 0);
    int32_t body = 0;
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
    scalar compressionRatio = 1.;

    vec emitterVelocity = vec(0,0);

    bool velocityNoise = false;
    bool densityNoise = false;
    bool areaNoise = false;
    scalar noiseAmplitude = 1.;
    scalar noiseFrequency = 1.;
    int32_t noiseOctaves = 4;
    int32_t noiseSeed = 0x12345678;

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
    bool inletOnce = true;
    std::size_t cellsX, cellsY;
    vec domainMin, domainMax;
    scalar baseSupport, baseRadius;


    YAML::Node config;
    void initializeSPH();
    void initializeParameters();

    std::vector<vec> genParticles(vec minCoord, vec maxCoord, scalar radius, scalar packing = packing_2D);
    std::vector<vec> genParticlesCircle(vec minCoord, vec maxCoord, scalar radius, scalar packing = packing_2D);
    template<typename Func> auto boundaryFunc(const int32_t ptclIndex, Func&& c) {
    auto p = fluidPosition[ptclIndex];
    auto h = fluidSupport[ptclIndex];
    for (auto ti : fluidTriangleNeighborList[ptclIndex]) {
        const auto& tri = boundaryTriangles[ti];
        auto [hit, pb, d, k, gk] = interactTriangle(p, h, tri);
        if (hit)
            c(pb, d, k, gk, true, ti);
    
}}
        template<typename Func> auto boundaryFunc(const vec& p, const scalar& h, Func&& c) {
            for(int32_t t = 0; t < boundaryTriangles.size(); ++t){
        auto tri = boundaryTriangles[t];
        auto [hit, pb, d, k, gk] = interactTriangle(p, h, tri);
        if (hit)
            c(pb, d, k, gk, true, t);
    }
        }


    std::vector<int32_t>& getCell(scalar x, scalar y);
    std::vector<int32_t>& getBoundaryCell(int32_t xi, int32_t yi);
    std::vector<int32_t>& getTriangleCell(int32_t xi, int32_t yi);

template<typename T> auto syncGhost(std::vector<T>& ptr){
    #pragma omp parallel for
    for(int32_t i = 0; i < ghostParticles.size(); ++i){
        int32_t idx = ghostParticles[i];
        ptr[idx] = ptr[fluidGhostIndex[idx]];
    }
}


public:
    SPHSimulation(std::string _config = ""){
        if(_config != ""){
            config = YAML::Load(_config);
            initializeParameters();
            initializeSPH();
        }else{
            config = YAML::Load(
R"(
fluids:
    - min: [0.2, 0.4]
      max: [0.4, 0.6]
      velocity: [1,0]
      radius: 0.00279783
      type: once
      shape: spherical
    - min: [0.6, 0.4]
      max: [0.8, 0.6]
      velocity: [-1,0]
      radius: 0.00279783
      type: once
      shape: spherical

gravity:
    - pointSource: false
      direction: [0., -1.]
      magnitude: 9.81

triangles:
    - v0: [0.375, 0.02]
      v1: [0.40625, 0.02]
      v2: [0.40625, 0.06]
    - v0: [0.40625, 0.02]
      v1: [0.59375, 0.02]
      v2: [0.59375, 0.06]
    - v0: [0.40625, 0.02]
      v1: [0.40625, 0.06]
      v2: [0.59375, 0.06]
    - v0: [0.59375, 0.02]
      v1: [0.625, 0.02]
      v2: [0.59375, 0.06]

domain:
    min: [0,0]
    max: [1, 1]
    epsilon: 0.02

sim:
    maxDt: 0.0005
    minDt: 0.0005

export:
    active: true
    limit: 2.0
    interval: 5
)"

);
            initializeParameters();
            initializeSPH();
        }
    }

    std::vector<vec> fluidPosition, fluidVelocity, fluidAccel, fluidPredVelocity, fluidPredAccel, fluidPredPosition;
    std::vector<scalar> fluidDensity, fluidVorticity, fluidAngularVelocity, fluidArea, fluidRestDensity, fluidSupport, fluidInitialAngularVelocity;
    std::vector<scalar> fluidAlpha, fluidActualArea, fluidPressure1, fluidPressure2, fluidBoundaryPressure, fluidSourceTerm, fluidDpDt, fluidDensityStar, fluidPriorPressure;
    std::vector<std::vector<int32_t>> fluidNeighborList, fluidTriangleNeighborList;
    std::vector<int32_t> fluidUID;
    std::vector<int32_t> fluidGhostIndex;
    int32_t fluidCounter = 0;
    std::vector<std::tuple<scalar, scalar, scalar>> boundaryBarycentricPressure;
    std::vector<vec> boundaryPolygon;
    std::vector<std::vector<vec>> boundaryObstacles;

    std::vector<std::vector<int32_t>> cellArray, cellBoundaryArray, cellTriangleArray;

    std::vector<std::pair<vec,vec>> boundaryParticles;

    std::vector<Triangle> boundaryTriangles;
    std::vector<boost::atomic<scalar>*> boundaryPressureForceX, boundaryPressureForceY, boundaryDragForceX, boundaryDragForceY;

    std::vector<fluidSource> fluidSources;
    std::vector<gravitySource> gravitySources;
    std::vector<int32_t> ghostParticles;

    std::vector<vec> fluidColorFieldGradient, fluidColorFieldGradientDifference, fluidColorFieldGradientSymmetric;
    std::vector<scalar> fluidColorField;

    std::vector<vec> fluidAdvectionVelocity, fluidPressureAccel, fluidPressureAccelSimple, fluidPressureAccelDifference, fluidPressureVelocity, fluidInitialPosition, fluidInitialVelocity;

    enum struct property_t{
        position, velocity, accel, predVelocity, predAccel, predPosition,
        density, vorticity, angularVelocity, area, restDensity, support,
        alpha, actualArea, pressure1, pressure2, boundaryPressure, sourceTerm, dpdt, rhoStar, priorPressure,
        UID, neighbors, ghostIndex,
        color, colorGrad, colorGradSymm, colorGradDiff
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

    void evalPressureForce();

    void XSPH();
    void computeVorticity();
    void refineVorticity();
    void externalForces();
    void BXSPH();

    void gaseousSPH();

    void dump();

    std::tuple<scalar, scalar, std::vector<scalar>> colorMap(property_t prop, bool autoMinMax = true, scalar min = 0., scalar max = 1.);
   
    ParameterManager pm;

    fs::path resolveFile(std::string fileName, std::vector<std::string> search_paths = {});

    scalar getScalar(std::string property){return pm.get<scalar>(property);}
    vec getVec(std::string property){return pm.get<vec>(property);}
    int32_t getInteger(std::string property){return pm.get<int32_t>(property);}
    bool getBoolean(std::string property){return pm.get<bool>(property);}

    void setScalar(std::string property, scalar value){pm.get<scalar>(property) = value;}
    void setVec(std::string property, vec value){pm.get<vec>(property) = value;}
    void setInteger(std::string property, int32_t value){pm.get<int32_t>(property) = value;}
    void setBoolean(std::string property, bool value){pm.get<bool>(property) = value;}


};