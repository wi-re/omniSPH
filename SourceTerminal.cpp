#include <tools/exceptionHandling.h>
#include <sstream>
#include <thread>
//#include <windows.h>
#include <iostream> 
#include <simulation/2DMath.h>
#include <iomanip>

#include <tools/timer.h>

int main(int argc, char* argv[]) 
//try 
{
//     simulationState = SPHSimulation(R"(
// fluids:
//     - min: [0.2, 0.2]
//       max: [0.4, 0.4]
//       radius: 0.0025
//       velocity: [1, 0]
//       shape: spherical
//       type: inlet
//       timeLimit: 1.0
//       ramp: 0.5
// fluidsB:
//     - min: [0.6, 0.6]
//       max: [0.8, 0.8]
//       radius: 0.00125
//       velocity: [-1, 0]
//       shape: spherical
//       type: once
    
// triangles:
//     - v0: [0.375, 0.02]
//       v1: [0.40625, 0.02]
//       v2: [0.40625, 0.06]
//     - v0: [0.40625, 0.02]
//       v1: [0.59375, 0.02]
//       v2: [0.59375, 0.06]
//     - v0: [0.40625, 0.02]
//       v1: [0.40625, 0.06]
//       v2: [0.59375, 0.06]
//     - v0: [0.59375, 0.02]
//       v1: [0.625, 0.02]
//       v2: [0.59375, 0.06]

// gravity:
//     - pointSource: true
//       location: [0.5, 0.75]
//       magnitude: 9.81
//     - pointSource: true
//       location: [0.5, 0.25]
//       magnitude: 9.81

// domain:
//     min: [0,0]
//     max: [2, 1]

// sim:
//     maxDt: 0.001
//     minDt: 0.001
//     )");

//     simulationState = SPHSimulation(R"(
// fluids:
//     - min: [0.0225, 0.0225]
//       max: [0.25, 0.5]
//       radius: 0.0025
//       type: once
//     - min: [0.75, 0.021]
//       max: [0.98, 0.5]
//       radius: 0.0025
//       type: once

// gravity:
//     - pointSource: false
//       direction: [0., -1.]
//       magnitude: 9.81

// triangles:
//     - v0: [0.375, 0.02]
//       v1: [0.40625, 0.02]
//       v2: [0.40625, 0.06]
//     - v0: [0.40625, 0.02]
//       v1: [0.59375, 0.02]
//       v2: [0.59375, 0.06]
//     - v0: [0.40625, 0.02]
//       v1: [0.40625, 0.06]
//       v2: [0.59375, 0.06]
//     - v0: [0.59375, 0.02]
//       v1: [0.625, 0.02]
//       v2: [0.59375, 0.06]

// domain:
//     min: [0,0]
//     max: [1, 1]
//     epsilon: 0.02

// sim:
//     maxDt: 0.001
//     minDt: 0.001

// export:
//     active: true
//     limit: 0.5
//     interval: 10
// )");

  if(argc == 1){
    throw std::invalid_argument("Not enough parameters");

  }
  std::string fileName = argv[1];
  std::ifstream t(fileName);
  std::stringstream buffer;
  buffer << t.rdbuf();  



   auto simulationState = SPHSimulation(buffer.str());
    simulationState.pm.get<bool>("export.active") = true;
    auto limit = simulationState.pm.get<scalar>("export.limit");
    if(limit < 0.)
    while(true)
      simulationState.timestep();
    
    while(simulationState.pm.get<scalar>("sim.time") <= limit)
      simulationState.timestep();


    for (auto tptr : TimerManager::getTimers()) {
        auto& t = *tptr;
        auto stats = t.getStats().value();
        scalar sum = 0.;
        for (auto s : t.getSamples()) {
            sum += s;
        }
        std::cout << std::setprecision(4);
        std::cout << t.getDecriptor() << ": " << stats.median << "ms / " << stats.avg << "ms @ " << stats.stddev << " : " << stats.min << "ms / " << stats.max << "ms, total: " << sum/1000. << "s\n";

    }

    return 0;
}
//CATCH_DEFAULT