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
auto fpc = R"(triangles:
    - v0: [ 0. ,  0. ]
      v1: [ -0.5 ,  0.25 ]
      v2: [ -0.5 ,  -0.25 ]
      body: 0
t2:
    - v0: [ 1. ,  -0.5 ]
      v1: [ 1.5 ,  -0.75 ]
      v2: [ 1.5 ,  -0.25 ]
      body: 1
    - v0: [ 1. ,  0.5 ]
      v1: [ 1.5 ,  0.75 ]
      v2: [ 1.5 ,  0.25 ]
      body: 2
fluids:
     - min: [-2, -1.]
       max: [4., 1.]
       radius: 0.01
       type: once
       velocity: [1.,0.]
       compression: 1.60011
       velocityNoise: false
       areaNoise: false
       noiseAmplitude: 2.0
       noiseOctaves: 1
       noiseFrequency: 1.0
       noiseSeed: 0
       density: 1000
     - min: [-2, -1.]
       max: [-1.9, 1.]
       radius: 0.01
       type: source
       velocity: [1.,0]
     - min: [3.9, -1.]
       max: [4.0, 1.]
       radius: 0.01
       type: source
       velocity: [1.,0]
  

       


gravity:
    - pointSource: true
      location: [.0, .0]
      magnitude: 0.0

triangless:
    - v0: [-1.5, 0.25]
      v1: [-1.0, 0.125]
      v2: [-1.0, 0.375]

video:
  active: false
  fps: 100.0

props:
    maxnumptcls: 512000
    backgroundPressure: true

dfsph:
  divergenceSolve: false
  densityEta: 0.001

domain:
    min: [-1, -1]
    max: [1, 1]
    epsilon: 0.02
    periodicX: true
    periodicY: true

sim:
    maxDt: 0.0025
    minDt: 0.0025
    incompressible: false
    kappa: 1.5

ptcl:
    boundaryViscosity: 0.5
    viscosityConstant: 0.01
vorticity:
    nu_t: 0.05
    angularViscositys: 0.01

export:
    active: false
    limit: 4.0
    interval: 1

colorMap:
  buffer: 'velocity'
  auto: true
  min: 0.
  max: 10.)";

auto noiseSimulation = R"(
fluids:
     - min: [-1, -1.]
       max: [1., 1.]
       radius: 0.01
       type: once
       velocity: [0.,0.]
       compression: 1.50005
       velocityNoise: true
       areaNoise: false
       noiseAmplitude: 1.0
       noiseOctaves: 8
       noiseFrequency: 2.0
       noiseSeed: 1337
       density: 1000
f:
     - min: [-0.1, -0.1]
       max: [0.1, 0.1]
       radius: 0.01
       type: source
       velocity: [1,0]
       compression: 1.
       ramp: 1.0


gravity:
    - pointSource: true
      location: [.0, .0]
      magnitude: 0.0

triangless:
    - v0: [-1.5, 0.25]
      v1: [-1.0, 0.125]
      v2: [-1.0, 0.375]

video:
  active: true
  fps: 100.0

props:
    maxnumptcls: 512000
    backgroundPressure: false

dfsph:
  divergenceSolve: true
  densityEta: 0.001

domain:
    min: [-1, -1]
    max: [1, 1]
    epsilon: 0.02
    periodicX: true
    periodicY: true

sim:
    maxDt: 0.0025
    minDt: 0.0025
    incompressible: false
    kappa: 1.5

ptcl:
    viscosityConstant: 0.05
vorticity:
    nu_t: 0.05

export:
    active: true
    limit: 2.5
    interval: 1

colorMap:
  buffer: 'velocity'
  auto: true
  min: 0.
  max: 10.
)";
    simulationState = SPHSimulation(fpc);
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