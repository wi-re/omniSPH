#pragma once
#define _USE_MATH_DEFINES
#include <array>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <tools/ParameterManager.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

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
constexpr inline scalar domainScale = 20.0;

inline scalar pointScale = 2.0 * domainScale;

// 1.000 compression
//constexpr inline scalar scale = 0.053599060968293831; // 100 rows
////constexpr inline scalar scale = 0.10607382169550233; // 50 rows
//constexpr inline scalar packing_2D = (scalar)0.21324955665222379 * scale;
//constexpr inline scalar spacing_2D = (scalar)0.19775018866158592 * scale;

// 0.9990
//constexpr inline scalar scale = 0.0.10606839266103468; // 50 rows
//constexpr inline scalar packing_2D = (scalar)0.21315991302770557 * scale;
//constexpr inline scalar spacing_2D = (scalar)0.19778427897675366 * scale;

constexpr inline scalar scale = 0.053599060968293831; // desired particles: 100
constexpr inline scalar packing_2D = 0.21314955665222379 * scale; // desired compression: 1
constexpr inline scalar spacing_2D = 0.19775018866158592 * scale; // actual delta: 0.99799999843898057

// 1.005 compression
//constexpr inline scalar scale = 0.10634501612361064; // 50 rows
//constexpr inline scalar packing_2D = (scalar)0.21263365441208518 * scale;
//constexpr inline scalar spacing_2D = (scalar)0.19604825172931667 * scale;

inline scalar inletPacking = packing_2D;
inline scalar inletSpacing = spacing_2D;



constexpr inline scalar offset = 0.125 / domainScale;
constexpr inline scalar minDt = 0.0001;
constexpr inline scalar maxDt = 0.002;
//constexpr inline scalar scale = 1.00;
constexpr inline scalar scale_2 = 1.4142135623730;
constexpr inline scalar pi = 3.141592653589793238;

// parameters for the overall simulation
// used to define the overall simulation domain (including the outer gap)

constexpr inline scalar domainWidth = (scalar)250.f / domainScale;
constexpr inline scalar domainHeight = (scalar)50.f / domainScale;

inline int32_t screenWidth = (int32_t) (domainWidth * 20.0 * domainScale);
inline int32_t screenHeight = (int32_t) (domainHeight * 20.0 * domainScale*2.0);
// used to create a gap from the domain to the edge of the window
constexpr static scalar domainEpsilon = (scalar)5.f / domainScale;
// used for direct forcing boundaries to reflect the velocity and reduce it's intensity
constexpr static scalar dampingFactor = (scalar)-0.5f;
// number of cells required for the cell grid. due to H = 1.0 this is just the domain size
constexpr inline std::size_t cellsX = std::size_t(domainWidth / scale);
constexpr inline std::size_t cellsY = std::size_t(domainHeight / scale);
// values indicating the actual domain size (excluding the outer gap)
inline vec domainMin(domainEpsilon, domainEpsilon);
inline vec domainMax(domainWidth - domainEpsilon, domainHeight - domainEpsilon);
// cell grid used for acceleration
inline std::array<std::vector<int32_t>, cellsX * cellsY> cellArray;

// parameters for particle properties
// fixed timestep
inline scalar dt = (scalar)0.0008;
// packing factor for a dense hexagonal grid of particles
// 0.21314955649168332
// 0.21294955649168332
// fixed particle radius such that H = 1.0
constexpr inline scalar radius = (scalar)0.22360679774997896 * scale; 
constexpr inline scalar support = (scalar) scale;
// area determined as pi radius^2
constexpr inline scalar area = pi * radius * radius;
// rest density of water
constexpr inline scalar rho0 = (scalar)998.0;
// helper value to simplify some terms
constexpr inline scalar mass = area * rho0;
// external gravity force
inline vec gravity((scalar)0.0, (scalar)-9.8);
//inline vec gravity((scalar)0.0, (scalar)0);
// factor used for XSPH
constexpr inline scalar viscosityConstant = (scalar)0.02;

// kernel Parameters
constexpr inline int32_t targetNeighbors = 20;
constexpr inline scalar kernelSize = (scalar)1.778002;
constexpr inline scalar kernelConstant = (scalar)(80.0 / (7.0 * M_PI * support * support));
constexpr inline scalar gradConstant = (scalar)(80.0 / (7.0 * M_PI * support * support * support));

// numerical Parameters
constexpr inline scalar epsilon((scalar)1e-7);

// Simple structure representing most particle quantities
// For performance a struct of arrays design would be better but this is not
// necessary for this simple simulaiton.
struct Particle {
    Particle(scalar _x, scalar _y) : pos(_x, _y), vel(0.f, 0.f), rho(0) { reset(); }
  inline void reset() {
    accel = vec(0, 0);
    rho = scalar(0.0);
    neighbors.clear();
  }
  vec pos = vec(0, 0);
  vec vel = vec(0, 0);
  vec accel = vec(0, 0);
  scalar rho = scalar(0);
  std::vector<int32_t> neighbors;
  scalar vorticity = scalar (0.0);
  scalar angularVelocity = scalar(0.0);
  int64_t uid = 0;
};
// Structure containing additional quantities required for a DFSPH-like pressure solver
struct dfsphState {
  inline void reset() {
    alpha = 0.0;
    area = 0.0;
    vel = vec(0, 0);
    pressure1 = 0.0;
    pressure2 = 0.0;
	pressureBoundary = 0.0;
    source = 0.0;
    accel = vec(0, 0);
    dpdt = rhoStar = 0.0;
  }
  scalar alpha = scalar(0);
  scalar area = scalar(0);
  scalar pressure1 = scalar(0);
  scalar pressure2 = scalar(0);
  scalar pressureBoundary = scalar(0);
  scalar source = scalar(0);
  scalar dpdt = scalar(0);
  scalar rhoStar = scalar(0);
  vec vel = vec(0, 0);
  vec accel = vec(0, 0);
};
struct Triangle {
  vec v0 = vec(0, 0), v1 = vec(0, 0), v2 = vec(0, 0);
};

// global simulation state for ease of use
inline std::vector<Particle> particles;
inline std::vector<dfsphState> particlesDFSPH;
inline std::vector<Triangle> triangles;
inline std::vector<vec> polygon;
inline  std::vector < std::vector<vec>> obstacles;
inline scalar simulationTime = 0.0;

// this function generates a hexagonal grid of particles from minCoord to maxCoord
// DOES NOT CHECK FOR PARTICLES OCCUPYING THIS SPACE
std::vector<Particle> genParticles(vec minCoord, vec maxCoord, scalar packing = packing_2D);
// initializes the SPH simulation with the given scene number
void initializeSPH(int32_t scene = 0);
// execute a single SPH timestep
void timestep();
void render();
void initRender();
// print all information about the given particle index for debugging
void printParticle(int32_t idx);
// converts a std::chrono duration to a millisecond value
scalar toMs(clk::duration dur);


void emitParticles();

enum struct cornerAngle {
	acute, ortho, obtuse, nobtuse, northo, nacute, Box4, Box1, Box1_4
};
enum struct boundaryMethod {
	analytical, semi, sdf
};

constexpr inline cornerAngle simulationCase = cornerAngle::Box4;
constexpr inline boundaryMethod simulationMethod = boundaryMethod::sdf;




std::vector<int32_t>& getCell(scalar x, scalar y);
std::pair<int32_t, int32_t> getCellIdx(scalar x, scalar y);
std::vector<int32_t>& getCell(int32_t xi, int32_t yi);