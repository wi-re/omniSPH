
#define _CRT_SECURE_NO_WARNINGS
#include "SPH.h"
#include "2DMath.h"
#include "config.h"
#include <algorithm>
#include <array>
#include <boost/range/combine.hpp>
#include <iostream>
#include <numeric>
#include <chrono>
#include <sstream>
#include <atomic>


#include <tools/timer.h>
#include <filesystem>

struct simulationData {
    std::byte* internalData = nullptr;

    // parameter based quantities
    std::size_t num_particles;
    int32_t frame, divergenceIterations, densityIterations;
    scalar time, timestep, divergenceError, densityError, radius;

    // computed data
    scalar minVelocity = DBL_MAX, maxVelocity = -DBL_MAX;
    scalar minAccel = DBL_MAX, maxAccel = -DBL_MAX;

    scalar minPressure = DBL_MAX, maxPressure = -DBL_MAX;
    scalar minDensity = DBL_MAX, maxDensity = -DBL_MAX;
    scalar minAngular = DBL_MAX, maxAngular = -DBL_MAX;
    scalar totalKineticEnergy = 0.0, totalPotentialEnergy = 0.0;

    // data
    vec* positions, * velocities, * accelerations;
    scalar* density, * pressure, * angularVelocity;
    scalar* kineticEnergy, * potentialEnergy;
    int64_t* UIDs;
};
#include <time.h>
#include <iomanip>
void dump() {
    auto& pm = ParameterManager::instance();
    if (pm.get<int32_t>("sim.frame") == 0) return;
    namespace fs = std::filesystem;
    static fs::path basePath = fs::current_path() / "data";
    static fs::path actualPath;

    static bool init = false;
    if (!init) {
        const std::chrono::time_point<std::chrono::system_clock> now =
            std::chrono::system_clock::now();
        const std::time_t t_c = std::chrono::system_clock::to_time_t(now);
        auto time = std::time(nullptr);
        std::stringstream ss;
        char buf[256];
        ss << std::put_time(std::localtime(&t_c), "%F_%T"); // ISO 8601 without timezone information.
        auto s = ss.str();
        std::replace(s.begin(), s.end(), ':', '-');
        actualPath = basePath / s;

        fs::create_directories(actualPath);
        auto summaryPath = actualPath / std::string("summary.sphlog");
        std::ofstream outFile;
        summaryFile.open(summaryPath, std::ios::out | std::ios::binary);
        summaryFileOpen = true;
        init = true;
    }
    std::size_t num_ptcls = particles.size();
    std::size_t payloadSize = sizeof(vec) * 3 + sizeof(scalar) * 5 + sizeof(int64_t);
    std::byte* rawData = (std::byte*)malloc(num_ptcls * payloadSize);

    simulationData data{ 
        .internalData = rawData, 
        .num_particles = num_ptcls,
        .frame = pm.get<int32_t>("sim.frame"), 
        .divergenceIterations = pm.get<int32_t>("dfsph.divergenceIterations"), .densityIterations = pm.get<int32_t>("dfsph.densityIterations"),
        .time = pm.get<scalar>("sim.time"), .timestep = pm.get<scalar>("sim.dt"),
        .divergenceError = pm.get<scalar>("dfsph.divergenceError"),.densityError = pm.get<scalar>("dfsph.densityError"),
        .radius = pm.get<scalar>("ptcl.radius") };

    std::cout << "Number of particles: " << data.num_particles << std::endl;
    std::cout << "Time: " << data.time << " ( dt = " << data.timestep << " ) @ frame " << data.frame << std::endl;
    std::cout << "Divergence/Density Iterations: " << data.divergenceIterations << " / " << data.densityIterations << std::endl;
    std::cout << "Divergence/Density Error: " << data.divergenceError << " / " << data.densityError << std::endl;

    for (int32_t i = 0; i < num_ptcls; ++i) {
        const auto& p = particles[i];
        const auto& dp = particlesDFSPH[i];

        auto update = [](auto& minEntry, auto& maxEntry, auto candidate) {
            minEntry = std::min(minEntry, candidate);
            maxEntry = std::max(maxEntry, candidate);
        };
        update(data.minVelocity, data.maxVelocity, p.vel.norm());
        update(data.minAccel, data.maxAccel, p.accel.norm());
        update(data.minPressure, data.maxPressure, dp.pressure1);
        update(data.minDensity, data.maxDensity, p.rho);
        update(data.minAngular, data.maxAngular, p.angularVelocity);
        data.totalKineticEnergy += 0.5 * mass * p.vel.squaredNorm();
        data.totalPotentialEnergy += gravitySwitch ? mass * gravity.norm() * (p.pos.y() - domainEpsilon) : 0.0;
    }
    std::cout << "Range of " << "velocity" << " : " << data.minVelocity << " -> " << data.maxVelocity << std::endl;
    std::cout << "Range of " << "accel" << " : " << data.minAccel << " -> " << data.maxAccel << std::endl;
    std::cout << "Range of " << "pressure" << " : " << data.minPressure << " -> " << data.maxPressure << std::endl;
    std::cout << "Range of " << "density" << " : " << data.minDensity << " -> " << data.maxDensity << std::endl;
    std::cout << "Range of " << "angular" << " : " << data.minAngular << " -> " << data.maxAngular << std::endl;
    std::cout << "Kinetic   Energy: " << data.totalKineticEnergy << std::endl;
    std::cout << "Potential Energy: " << data.totalPotentialEnergy << std::endl;
    std::cout << "Total Energy: " << data.totalKineticEnergy + data.totalPotentialEnergy << std::endl;

    std::cout << "rawData allocation size: " << payloadSize * num_ptcls << " : " << num_ptcls << " x " << payloadSize << std::endl;
    data.positions = (vec*)(rawData);
    data.velocities = (vec*)(rawData + sizeof(vec) * num_ptcls);
    data.accelerations = (vec*)(rawData + 2 * sizeof(vec) * num_ptcls);
    data.density = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls);
    data.pressure = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + sizeof(scalar) * num_ptcls);
    data.angularVelocity = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 2 * sizeof(scalar) * num_ptcls);
    data.kineticEnergy = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 3 * sizeof(scalar) * num_ptcls);
    data.potentialEnergy = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 4 * sizeof(scalar) * num_ptcls);
    data.UIDs = (int64_t*)(rawData + 3 * sizeof(vec) * num_ptcls + 5 * sizeof(scalar) * num_ptcls);



    for (int32_t i = 0; i < num_ptcls; ++i) {
        const auto& p = particles[i];
        const auto& dp = particlesDFSPH[i];

        data.positions[i] = p.pos;
        data.velocities[i] = p.vel;
        data.accelerations[i] = p.accel;
        data.density[i] = p.rho;
        data.pressure[i] = dp.pressure1;
        data.angularVelocity[i] = p.angularVelocity;
        data.kineticEnergy[i] = 0.5 * mass * p.vel.norm();
        data.potentialEnergy[i] = gravitySwitch ? mass * gravity.norm() * (p.pos.y() - domainEpsilon) : 0.0;
        data.UIDs[i] = p.uid;
    }

    std::stringstream sstream;
    sstream << std::setfill('0') << std::setw(5) << pm.get<int32_t>("sim.frame") << ".sph";

    fs::path file = actualPath / sstream.str();


    std::cout << "Base Path: " << basePath << std::endl;
    std::cout << "Actual Path: " << actualPath << std::endl;
    std::cout << "File: " << file << std::endl;

    std::ofstream outFile;
    outFile.open(file, std::ios::out | std::ios::binary);

    outFile.write(reinterpret_cast<char*>(&num_ptcls), sizeof(num_ptcls));
    outFile.write(reinterpret_cast<char*>(&data), sizeof(simulationData));
    outFile.write(reinterpret_cast<char*>(rawData), payloadSize * num_ptcls);

    summaryFile.write(reinterpret_cast<char*>(&data), sizeof(simulationData));

    outFile.close();
    free(rawData);
}


void timestep() {
    dt = 0.002;
    auto& speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    auto speedGoal = ParameterManager::instance().get<scalar>("sim.inletSpeedGoal");
    static auto& t = ParameterManager::instance().get<scalar>("sim.time");
    static bool once = true;
    if (once) {
        speed = 0.0;
        once = false;
    }
    if (t <= 1.0) {
        speed = speedGoal;
    }
    if (t > 7.0)
        inletSwitch = false;


    TIME_CODE(0, "Simulation - Overall",
        TIME_CODE(1, "Simulation - Reset", resetFrame());
        TIME_CODE(2, "Simulation - Cell construction", fillCells());
        TIME_CODE(3, "Simulation - Neighbor search", neighborList());
        TIME_CODE(4, "Simulation - Density", density());
        TIME_CODE(5, "Simulation - Vorticity", computeVorticity());
        TIME_CODE(6, "Simulation - External", externalForces());
        TIME_CODE(7, "Simulation - Divergence", divergenceSolve());
        TIME_CODE(8, "Simulation - Density", densitySolve());
        TIME_CODE(9, "Simulation - XSPH", XSPH());
        TIME_CODE(10, "Simulation - Vorticity", refineVorticity());
        TIME_CODE(11, "Simulation - Integration", Integrate());
        TIME_CODE(12, "Simulation - Dump", dump());
        TIME_CODE(13, "Simulation - Emission", emitParticles());
    );

    ParameterManager::instance().get<int32_t>("sim.frame")++;
  //auto cellBegin = clk::now();
  //fillCells();
  //auto neighBegin = clk::now();
  //neighborList();
  //auto densityBegin = clk::now();
  //density();
  //auto externalBegin = clk::now();
  //externalForces();
  //auto divSolveBegin = clk::now();
  //int32_t div_it = divergenceSolve();
  //auto incSolveBegin = clk::now();
  //int32_t den_it = densitySolve();
  //auto xsphBegin = clk::now();
  //XSPH();
  //auto integrateBegin = clk::now();
  //Integrate();
  //auto endFrame = clk::now();
  //std::stringstream sstream;
  //sstream << "########################################################################################################\n";
  ////sstream << "Frame took " << toMs(endFrame - beginFrame) << "ms particle count " << particles.size() << "\n";
  //sstream << "Cells      " << toMs(neighBegin - cellBegin) << "ms Neighbors      " << toMs(densityBegin - neighBegin) << "\n";
  //sstream << "Density    " << toMs(externalBegin - densityBegin) << "ms External       " << toMs(divSolveBegin - externalBegin) << "\n";
  //sstream << "Divergence " << toMs(incSolveBegin - divSolveBegin) << "ms Incompressible " << toMs(xsphBegin - incSolveBegin) << "\n";
  //sstream << "XSPH       " << toMs(integrateBegin - xsphBegin) << "ms Integrate      " << toMs(endFrame - integrateBegin) << "\n";
  //sstream << "Divergence Iterations " << div_it << ", Incompressibility Iterations " << den_it << "\n";
  //std::cout << sstream.str();
}

using clk = std::chrono::high_resolution_clock;
scalar toMs(clk::duration dur) {
  return static_cast<scalar>(std::chrono::duration_cast<std::chrono::microseconds>(dur).count()) / scalar(1000.0);
}

std::vector<Particle> genParticles(vec minCoord, vec maxCoord, scalar packing) {
  auto gen_position = [](auto r, auto i, auto j) -> vec {
    return vec(r * (2.0 * scalar(i) + scalar(j % 2)), r * ::sqrt(3.0) * scalar(j));
  };
  auto diff = maxCoord - minCoord;
  auto requiredSlices_x = ::ceil(diff.x() / packing);
  auto requiredSlices_y = ::ceil(diff.y() / (::sqrt(3.0) * packing));
 // std::cout << "Generating particles on a " << requiredSlices_x << " x " << requiredSlices_y << " hex grid " << std::endl;
  std::vector<Particle> points;
  for (int32_t x_it = 0; x_it < requiredSlices_x + 1; ++x_it)
    for (int32_t y_it = 0; y_it < requiredSlices_y + 1; ++y_it) {
      vec p = minCoord;
      vec g = gen_position(packing, x_it, y_it);
      vec pos = p + g;
      if (pos.x() < maxCoord.x() && pos.y() <= maxCoord.y()+ scale * 0.2)
        points.emplace_back(pos.x(), pos.y());
    }
  return points;
}

void printParticle(int32_t idx) {
  auto &p = particles[idx];
  auto &d = particlesDFSPH[idx];
  std::cout << "###############################################\n";
  std::cout << "Particle: " << idx << "\n";
  std::cout << "rho: " << p.rho << "\n";
  std::cout << "pos: [" << p.pos.x() << " : " << p.pos.y() << "], vel: [" << p.vel.x() << " : " << p.vel.y()
            << "], acc: [" << p.accel.x() << " : " << p.accel.y() << "]\n";
  std::cout << "alpha: " << d.alpha << ", area: " << d.area << ", source: " << d.source << ", boundaryP: " << d.pressureBoundary <<"\n";
  std::cout << "pressure1: " << d.pressure1 << ", pressure2: " << d.pressure2 << ", dpdt: " << d.dpdt
            << ", rho*: " << d.rhoStar << "\n";
  std::cout << "predictedVel: [" << d.vel.x() << " : " << d.vel.y() << "], pressureAccel: [" << d.accel.x() << " : "
            << d.accel.y() << "]" << std::endl;
}
