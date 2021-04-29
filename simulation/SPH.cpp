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
        TIME_CODE(12, "Simulation - Emission", emitParticles());
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
