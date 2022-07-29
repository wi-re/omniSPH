#include <tools/exceptionHandling.h>
#include "gui/glui.h"
#include <sstream>
#include <thread>
//#include <windows.h>
#include <iostream> 
#include <simulation/2DMath.h>
#include <iomanip>

#include <tools/timer.h>

scalar sdpolygon2(std::vector<vec> v, vec p) {
    scalar d = (p - v[0]).dot(p - v[0]);
    scalar s = 1.0;
    for (int i = 0, j = v.size() - 1; i < v.size(); j = i, i++) {
        vec e = v[j] - v[i];
        vec w = p - v[i];
        vec b = w - e * std::clamp(w.dot(e) / e.dot(e), 0.0, 1.0);
        d = std::min(d, b.dot(b));
        bool b1 = p.y() >= v[i].y();
        bool b2 = p.y() < v[j].y();
        bool b3 = e.x() * w.y() > e.y() * w.x();
        if ((b1 && b2 && b3) || (!b1 && !b2 && !b3))
            s *= -1.0;
    }
    return s * sqrt(d);
}

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

// auto n = 5;
// std::vector<vec> positions;
// auto area = 1.;
// auto radius = sqrt(area / double_pi);
// auto support = sqrt(area * targetNeighbors / double_pi);
// auto packing = 0.39900743165053487;

// for(int32_t x = -n; x <= n; ++x){
// for(int32_t y = -n; y <= n; ++y){
//   positions.push_back(vec{
//     packing * support * (double) x,
//     packing * support * (double) y
//   });
// }}

// vec center{0.,0.};
// auto rho = 0.;
// for(auto p : positions)
//   rho += area * W(center, p, support, support);
// printf("Source.cpp: radius: %g, area: %g, support: %g, packing: %g\n", radius, area, support, packing * support);
// std::cout << "Actual density: " << rho << std::endl;
// fluidSource source{
//   .emitterMin = vec(-packing * support * 5.,-packing * support * 5.),
//   .emitterMax = vec(packing * support * 5.,packing* support  * 5.),
//   .emitterRadius = radius
// };
// auto ptcls = source.genParticles();
// rho = 0.;
// for(auto p : ptcls)
//   rho += area * W(center, p, support, support);
// std::cout << "Gen density: " << rho << std::endl;
// auto rSum = 0.;
// auto lSum = 0.;
// auto dSum = 0.;
// for(int32_t i = 0; i < positions.size(); ++i){
//   auto p = positions[i];
//   auto b = ptcls[i];
//   auto d = b - center;
//   printf("[%d]: [%g, %g] - [%g, %g]\n", i, p.x(), p.y(), b.x(), b.y());
//   printf("\t[%g, %g] -> %g -> Kernel (m): %g\n", (p-center).x(), (p-center).y(), (p-center).norm(), W(p, center, support, support));
//   printf("\t[%g, %g] -> %g -> Kernel (a): %g\n", (b-center).x(), (b-center).y(), (b-center).norm(), W(b, center, support, support));
//   printf("\tDiff: %g\n", W(p, center, support, support) - W(b, center, support, support));
//   rSum += W(p, center, support, support);
//   lSum += W(b, center, support, support);
//   dSum += abs(W(p, center, support, support) - W(b, center, support, support));
// }
// printf("%g : %g -> %g || %g\n", rSum, lSum, rSum - lSum, dSum);


if(argc > 1){
  std::string fileName = argv[1];
  std::ifstream t(fileName);
  std::stringstream buffer;
  buffer << t.rdbuf();  
    simulationState = SPHSimulation(buffer.str());
}
else{



    simulationState = SPHSimulation(R"(
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
)");
}
    auto& gui = GUI::instance();
    gui.render_lock.lock();
    gui.initSimulation();
    gui.initParameters(argc, argv);
    gui.initGVDB();
    gui.initGL(argc, argv);
    gui.renderLoop();


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