#include <simulation/SPH.h>

#include <cfloat>
#include <tools/timer.h>

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
    scalar minError = DBL_MAX, maxError = -DBL_MAX;
    scalar totalKineticEnergy = 0.0, totalPotentialEnergy = 0.0;
    scalar minKineticEnergy = DBL_MAX, maxKineticEnergy = -DBL_MAX;
    scalar minPotentialKineticEnergy = DBL_MAX, maxPotentialEnergy = -DBL_MAX;
	scalar timer_overall, timer_reset, timer_fill, timer_neighbors, timer_density, timer_computeVorticity, timer_external, timer_divergenceSolve, timer_densitySolve, timer_xsph, timer_refineVorticity, timer_integrate, timer_BXSPH, timer_dump, timer_emit;
};

void SPHSimulation::dump(){
//    auto& pm = pm;
    auto& interval = pm.get<int32_t>("export.interval");
if(lastExport == 0) lastExport = -interval;
    auto& t = pm.get<scalar>("sim.time");
    auto& f = pm.get<int32_t>("sim.frame");
    
    if(!pm.get<bool>("export.active")) return;

   if (pm.get<int32_t>("sim.frame") == 0) return;
   fs::path basePath = fs::current_path() / "data";


    if(!(f >= lastExport + interval)) return;

    lastExport = f;

   if (!initDump) {
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
       initDump = true;


       auto configPath = actualPath / std::string("config.triangles");
       outFile.open(configPath, std::ios::out);
       outFile << "Triangle count: " << boundaryTriangles.size() << std::endl;
       int32_t i = 0; 
       for (auto& t : boundaryTriangles) {
           outFile << i++ << " : " << t.v0.x() << " " << t.v0.y() << " - " << t.v1.x() << " " << t.v1.y() << " - " << t.v2.x() << " " << t.v2.y() << std::endl;
       }
       outFile.close();


      {
        auto configPath = actualPath / std::string("config.yaml");
        auto config = pm.buildTree();
       outFile.open(configPath, std::ios::out);
       outFile << config << std::endl;
       outFile.close();
}

       std::ofstream boundaryFile;
       auto boundaryPath = actualPath / std::string("boundary.particles");
       boundaryFile.open(boundaryPath, std::ios::out | std::ios::binary);
       // std::cout << 
       for(int32_t b = 0; b < boundaryParticles.size(); ++b){
           boundaryFile.write(reinterpret_cast<char*>(&boundaryParticles[b].first), sizeof(vec));
       }
       for(int32_t b = 0; b < boundaryParticles.size(); ++b){
           boundaryFile.write(reinterpret_cast<char*>(&boundaryParticles[b].second), sizeof(vec));
       }

       boundaryFile.close();
   }
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
   std::size_t num_ptcls = numPtcls;
//    std::size_t payloadSize = sizeof(vec) * 4 + sizeof(scalar) * 5 + sizeof(int64_t);
//    std::byte* rawData = (std::byte*)malloc(num_ptcls * payloadSize);

   simulationData data{ 
    //    .internalData = rawData, 
       .num_particles = num_ptcls,
       .frame = pm.get<int32_t>("sim.frame"), 
       .divergenceIterations = pm.get<int32_t>("dfsph.divergenceIterations"), .densityIterations = pm.get<int32_t>("dfsph.densityIterations"),
       .time = pm.get<scalar>("sim.time"), .timestep = pm.get<scalar>("sim.dt"),
       .divergenceError = pm.get<scalar>("dfsph.divergenceError"),.densityError = pm.get<scalar>("dfsph.densityError"),
       .radius = pm.get<scalar>("props.baseRadius") };

   std::cout << "ptcls: " << data.num_particles << "\t";
   std::cout << "t: " << data.time << " ( dt = " << data.timestep << " )\tframe " << data.frame << " \t";
   std::cout << "div/inc it: " << data.divergenceIterations << "/" << data.densityIterations << " \t";
   std::cout << "div/inc err: " << data.divergenceError << "/" << data.densityError << "\t->\t";

   for (int32_t i = 0; i < num_ptcls; ++i) {

       auto update = [](auto& minEntry, auto& maxEntry, auto candidate) {
           minEntry = std::min(minEntry, candidate);
           maxEntry = std::max(maxEntry, candidate);
       };
       update(data.minVelocity, data.maxVelocity, fluidVelocity[i].norm());
       update(data.minAccel, data.maxAccel, fluidAccel[i].norm());
       update(data.minPressure, data.maxPressure, fluidPressure1[i]);
       update(data.minDensity, data.maxDensity, fluidDensity[i]);
       update(data.minAngular, data.maxAngular, fluidAngularVelocity[i]);
       update(data.minError, data.maxError, fluidDensityStar[i] / fluidArea[i]);
       update(data.minKineticEnergy, data.maxKineticEnergy, 0.5 * fluidArea[i] * fluidRestDensity[i] * fluidVelocity[i].squaredNorm());

       update(data.minPotentialKineticEnergy, data.maxPotentialEnergy, fluidArea[i] * fluidRestDensity[i] * (-9.81) * (fluidPosition[i].y() - domainMin.y()));
       data.totalKineticEnergy += 0.5 * fluidArea[i] * fluidRestDensity[i] * fluidVelocity[i].squaredNorm();
       data.totalPotentialEnergy += fluidArea[i] * fluidRestDensity[i] * (-9.81) * (fluidPosition[i].y() - domainMin.y());
   }
  //  std::cout << "Range of " << "velocity" << " : " << data.minVelocity << " -> " << data.maxVelocity << std::endl;
  //  std::cout << "Range of " << "accel" << " : " << data.minAccel << " -> " << data.maxAccel << std::endl;
  //  std::cout << "Range of " << "pressure" << " : " << data.minPressure << " -> " << data.maxPressure << std::endl;
  //  std::cout << "Range of " << "density" << " : " << data.minDensity << " -> " << data.maxDensity << std::endl;
  //  std::cout << "Range of " << "angular" << " : " << data.minAngular << " -> " << data.maxAngular << std::endl;
  //  std::cout << "Kinetic   Energy: " << data.totalKineticEnergy << std::endl;
  //  std::cout << "Potential Energy: " << data.totalPotentialEnergy << std::endl;
  //  std::cout << "Total Energy: " << data.totalKineticEnergy + data.totalPotentialEnergy << std::endl;

//    std::cout << "rawData allocation size: " << payloadSize * num_ptcls << " : " << num_ptcls << " x " << payloadSize << std::endl;
//    data.positions = (vec*)(rawData);
//    data.velocities = (vec*)(rawData + sizeof(vec) * num_ptcls);
//    data.accelerations = (vec*)(rawData + 2 * sizeof(vec) * num_ptcls);
//    data.density = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls);
//    data.pressure = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + sizeof(scalar) * num_ptcls);
//    data.angularVelocity = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 2 * sizeof(scalar) * num_ptcls);
//    data.kineticEnergy = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 3 * sizeof(scalar) * num_ptcls);
//    data.potentialEnergy = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 4 * sizeof(scalar) * num_ptcls);
//    data.UIDs = (int64_t*)(rawData + 3 * sizeof(vec) * num_ptcls + 5 * sizeof(scalar) * num_ptcls);
//    data.error = (scalar*)(rawData + 3 * sizeof(vec) * num_ptcls + 5 * sizeof(scalar) * num_ptcls + sizeof(int64_t)*num_ptcls);


//    std::vector<std::pair<int64_t,int64_t>> indices;
//    for (int32_t i = 0; i < num_ptcls; ++i) {
//        const auto& p = particles[i];
//        indices.push_back(std::make_pair(i, p.uid));
//    }
//    std::sort(indices.begin(), indices.end(),[](const auto& lhs, const auto& rhs){ return lhs.second < rhs.second;});



//    for (int32_t i = 0; i < num_ptcls; ++i) {
//        auto idx = indices[i].first;
//        const auto& p = particles[idx];
//        const auto& dp = particlesDFSPH[idx];
//        //std::cout << i << " " << indices[i].first << " - " << indices[i].second << "\n";


//        data.positions[i] = p.pos;
//        data.velocities[i] = p.vel;
//        data.accelerations[i] = p.accel;
//        data.density[i] = p.rho;
//        data.pressure[i] = dp.pressure1;
//        data.angularVelocity[i] = p.angularVelocity;
//        data.kineticEnergy[i] = 0.5 * mass * p.vel.norm();
//        data.potentialEnergy[i] = gravitySwitch ? mass * gravity.norm() * (p.pos.y() - domainEpsilon) : 0.0;
//        data.UIDs[i] = p.uid;
//        data.error[i] = dp.rhoStar / area;
//    }
   auto& timers = TimerManager::getTimers();
   for (auto t : timers) {
       if (t->getDecriptor() == "Simulation - Overall")
           data.timer_overall = t->getSamples().size() != 0 ? *(t->getSamples().end()-1) : 0.;
       if (t->getDecriptor() == "Simulation - Reset")
           data.timer_reset = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Cell construction")
           data.timer_fill = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Neighbor search")
           data.timer_neighbors = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Density")
           data.timer_density = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Vorticity")
           data.timer_computeVorticity = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - External")
           data.timer_external = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Divergence")
           data.timer_divergenceSolve = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Density")
           data.timer_densitySolve = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - XSPH")
           data.timer_xsph = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Vorticity")
           data.timer_refineVorticity = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Integration")
           data.timer_integrate = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - BXSPH")
           data.timer_BXSPH = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Dump")
           data.timer_dump = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
       if (t->getDecriptor() == "Simulation - Emission")
           data.timer_emit = t->getSamples().size() != 0 ? *(t->getSamples().end() - 1) : 0.;
   }
   //std::cout << "Timers: " << std::endl;
   //std::cout << " " << data.timer_overall << std::endl;
   //std::cout << " " << data.timer_reset << std::endl;
   //std::cout << " " << data.timer_fill << std::endl;
   //std::cout << " " << data.timer_neighbors << std::endl;
   //std::cout << " " << data.timer_computeVorticity << std::endl;
   //std::cout << " " << data.timer_external << std::endl;
   //std::cout << " " << data.timer_divergenceSolve << std::endl;
   //std::cout << " " << data.timer_densitySolve << std::endl;
   //std::cout << " " << data.timer_xsph << std::endl;
   //std::cout << " " << data.timer_refineVorticity << std::endl;
   //std::cout << " " << data.timer_integrate << std::endl;
   //std::cout << " " << data.timer_BXSPH << std::endl;
   //std::cout << " " << data.timer_dump << std::endl;
   //std::cout << " " << data.timer_emit << std::endl;


   auto frame = exportFrameCounter++;
   std::stringstream sstream;
   sstream << std::setfill('0') << std::setw(5) << frame;
   // sstream << std::setfill('0') << std::setw(5) << pm.get<int32_t>("sim.frame") << ".sph";

   fs::path file = actualPath / sstream.str();
   file.replace_extension("sph");

  //  std::cout << "Base Path: " << basePath << std::endl;
  std::cout << "Path: " << file << std::endl;
  //  std::cout << "File: " << file << std::endl;

   std::ofstream outFile;
   outFile.open(file, std::ios::out | std::ios::binary);

   outFile.write(reinterpret_cast<char*>(&num_ptcls), sizeof(num_ptcls));
   outFile.write(reinterpret_cast<char*>(&data), sizeof(simulationData));
//    outFile.write(reinterpret_cast<char*>(rawData), payloadSize * num_ptcls);
    fs::create_directories(actualPath/sstream.str());
    auto writeData = [&](auto data, std::string extension){
        using T = std::decay_t<decltype(data[0])>;

        fs::path file = actualPath / sstream.str() / "data";
        file.replace_extension(extension);
        std::ofstream outFile;
        outFile.open(file, std::ios::out | std::ios::binary);
        outFile.write(reinterpret_cast<char*>(data.data()), sizeof(T) * numPtcls);
        outFile.close();
    };
    writeData(fluidInitialPosition, "position");
    writeData(fluidPosition, "finalPosition");
    writeData(fluidInitialVelocity, "velocity");
    writeData(fluidVelocity, "finalVelocity");
    // writeData(fluidAccel, "accel");
    writeData(fluidPressure1, "pressure");
    writeData(fluidDensity, "density");
    // writeData(fluidVorticity, "vorticity");
    writeData(fluidAngularVelocity, "finalAngularVelocity");
    writeData(fluidInitialAngularVelocity, "angularVelocity");
    writeData(fluidArea, "area");
    // writeData(fluidRestDensity, "restDensity");
    // writeData(fluidSupport, "support");
    writeData(fluidUID, "uid");

    writeData(fluidGhostIndex, "ghostIndex");
    // writeData(fluidAdvectionVelocity, "advectionVelocity");
    // writeData(fluidPressureAccel, "pressureAccelSymmetric");    
    // writeData(fluidPressureAccelSimple, "pressureAccel");  
    // writeData(fluidPressureAccelDifference, "pressureAccelDifference");  
    // writeData(fluidPressureVelocity, "pressureVelocity");

    writeData(fluidColorField, "colorField");    
    writeData(fluidColorFieldGradient, "colorFieldGradient");  
    // writeData(fluidColorFieldGradientDifference, "colorFieldGradientDifference");  
    // writeData(fluidColorFieldGradientSymmetric, "colorFieldGradientSymmetric");


   summaryFile.write(reinterpret_cast<char*>(&data), sizeof(simulationData));

   outFile.close();
//    free(rawData);

}



void SPHSimulation::resetFrame(){
    bool periodicX = pm.get<bool>("domain.periodicX");
    bool periodicY = pm.get<bool>("domain.periodicY");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    for(int32_t i = 0; i < numPtcls; ++i){
    fluidAccel[i] = vec(0.,0.);
    fluidPredVelocity[i] = vec(0.,0.);
    fluidPredAccel[i] = vec(0.,0.);
    fluidPredPosition[i] = vec(0.,0.);
    fluidDensity[i] = 0.;
    //fluidVorticity[i] = 0.;
    //fluidAngularVelocity[i] = 0.;
    fluidAlpha[i] = 0.;
    fluidActualArea[i] = 0.;
    fluidPressure1[i] = 0.;
    fluidPressure2[i] = 0.;
    fluidBoundaryPressure[i] = 0.;
    fluidSourceTerm[i] = 0.;
    fluidDpDt[i] = 0.;
    fluidDensityStar[i] = 0.;
    fluidNeighborList[i] = std::vector<int32_t>{};
    fluidTriangleNeighborList[i] = std::vector<int32_t>{};
    fluidGhostIndex[i] = -1;

    fluidPressureAccel[i]=vec(0,0);
    fluidPressureAccelDifference[i]=vec(0,0);
    fluidPressureAccelSimple[i]=vec(0,0);
    fluidPressureVelocity[i]=vec(0,0);
    fluidAdvectionVelocity[i]=vec(0,0);
    fluidInitialVelocity[i]=vec(0,0);

}
ghostParticles.clear();
for(int32_t c = 0; c < cellsX * cellsY; ++c){
    cellArray[c] = std::vector<int32_t>{};
}
for(int32_t i = 0; i < boundaryPressureForceX.size(); ++ i){
  *boundaryPressureForceX[i] = 0.;
  *boundaryPressureForceY[i] = 0.;
  *boundaryDragForceX[i] = 0.;
  *boundaryDragForceY[i] = 0.;
}
auto virtualMin = pm.get<vec>("domain.virtualMin");
auto virtualMax = pm.get<vec>("domain.virtualMax");

        auto domainEpsilon = pm.get<scalar>("domain.epsilon");
        std::vector<int32_t> filteredParticles;
        filteredParticles.reserve(numPtcls);
        auto& t = pm.get<scalar>("sim.time");
        int32_t filteredCount = 0;
        for (int32_t i = 0; i < numPtcls; ++ i) {
            bool filtered = false;
            for(const auto& source : fluidSources){
                if(source.emitter != emitter_t::outlet) continue;
                if (fluidPosition[i].x() > source.emitterMin.x() && fluidPosition[i].x() < source.emitterMax.x()&&
                fluidPosition[i].y() > source.emitterMin.y() && fluidPosition[i].y() < source.emitterMax.y()) {
                ++filteredCount;
                filtered = true;
                break;
            }

            }
            if (filtered)
              continue;
            auto [ci, cj] =
                getCellIdx(fluidPosition[i].x(), fluidPosition[i].y());

            if (!periodicX && !periodicY) {
              if (fluidPosition[i].x() < domainMin.x() + 0.9 * domainEpsilon ||
                  fluidPosition[i].x() > domainMax.x() - 0.9 * domainEpsilon ||
                  fluidPosition[i].y() < domainMin.y() + 0.9 * domainEpsilon ||
                  fluidPosition[i].y() > domainMax.y() - 0.9 * domainEpsilon) {
                ++filteredCount;
              } else {
                filteredParticles.push_back(i);
              }
            } else if (periodicX && !periodicY) {
              if (fluidPosition[i].x() < virtualMin.x() || fluidPosition[i].x() > virtualMax.x() ||
                  fluidPosition[i].y() < domainMin.y() + 0.9 * domainEpsilon ||
                  fluidPosition[i].y() > domainMax.y() - 0.9 * domainEpsilon) {
                ++filteredCount;
              } else {
                filteredParticles.push_back(i);
              }
            } else if (periodicX && !periodicX) {
              if (fluidPosition[i].x() < domainMin.x() + 0.9 * domainEpsilon ||
                  fluidPosition[i].x() > domainMax.x() - 0.9 * domainEpsilon ||
                  fluidPosition[i].y() < virtualMin.y() || fluidPosition[i].y() > virtualMax.y()) {
                ++filteredCount;
              } else {
                filteredParticles.push_back(i);
              }
            } else {
              if (fluidPosition[i].x() < virtualMin.x() || fluidPosition[i].x() > virtualMax.x() || 
                  fluidPosition[i].y() < virtualMin.y() || fluidPosition[i].y() > virtualMax.y()) {
                ++filteredCount;
              } else {
                filteredParticles.push_back(i);
              }
            }
        }
        if (filteredCount > 0) {
            numPtcls = numPtcls - filteredCount;
            for (int32_t i = 0; i < numPtcls; ++i) {
                auto srcIdx = filteredParticles[i];
                fluidPosition[i] = fluidPosition[srcIdx];
                fluidVelocity[i] = fluidVelocity[srcIdx];
                fluidAngularVelocity[i] = fluidAngularVelocity[srcIdx];
                fluidVorticity[i] = fluidVorticity[srcIdx];
                fluidPriorPressure[i] = fluidPriorPressure[srcIdx];
                fluidArea[i] = fluidArea[srcIdx];
                fluidRestDensity[i] = fluidRestDensity[srcIdx];
                fluidUID[i]=fluidUID[srcIdx];
              fluidInitialPosition[i] = fluidPosition[srcIdx];
                fluidSupport[i] = fluidSupport[srcIdx];
            }

 }
// create ghost particles
/*
The domain:
*-*-------------*-*
|A|      B      |C|
*-*-------------*-*
| |             | |
| |             | |
|H|      0      |D|
| |             | |
| |             | |
*-*-------------*-*
|G|      F      |E|
*-*-------------*-*
Is mapped as
*/
// fillCells();
// static int32_t oldGhostCount = 0;
int32_t ghostCount = 0;
auto emitGhost = [&numPtcls, &ghostCount, this](int32_t srcIdx, scalar xOffset, scalar yOffset){  
      fluidPosition[numPtcls] = fluidPosition[srcIdx] + vec(xOffset, yOffset);
      fluidVelocity[numPtcls] = fluidVelocity[srcIdx];
      fluidAngularVelocity[numPtcls] = fluidAngularVelocity[srcIdx];
      fluidVorticity[numPtcls] = fluidVorticity[srcIdx];
      fluidPriorPressure[numPtcls] = fluidPriorPressure[srcIdx];
      fluidArea[numPtcls] = fluidArea[srcIdx];
      fluidRestDensity[numPtcls] = fluidRestDensity[srcIdx];
      fluidUID[numPtcls]=fluidUID[srcIdx];
      fluidInitialPosition[numPtcls] = fluidPosition[numPtcls];
      fluidSupport[numPtcls] = fluidSupport[srcIdx];
      fluidGhostIndex[numPtcls] = srcIdx;
      ghostParticles.push_back(numPtcls);
      numPtcls++;
      ghostCount++;
};

auto fringe = baseSupport * (double) pm.get<int32_t>("domain.buffer");
if(!periodicX && !periodicY) return;
int32_t actualPtcls = numPtcls;
for(int32_t i = 0; i < actualPtcls; ++i){
  auto pi = fluidPosition[i];
  if(periodicX){
    if(pi.x() < virtualMin.x() + fringe){
      emitGhost(i, virtualMax.x() - virtualMin.x(), 0.);
    }
    if(pi.x() > virtualMax.x() - fringe){
      emitGhost(i, virtualMin.x() - virtualMax.x(), 0.);
    }
  }
  if(periodicY){
    if(pi.y() < virtualMin.y() + fringe){
      emitGhost(i,0., virtualMax.y() - virtualMin.y());
    }
    if(pi.y() > virtualMax.y() - fringe){
      emitGhost(i, 0.,virtualMin.y() - virtualMax.y());
    }
  }
  if(periodicX && periodicY){
    if(pi.x() < virtualMin.x() + fringe && pi.y() < virtualMin.y() + fringe)
      emitGhost(i, virtualMax.x() - virtualMin.x(), virtualMax.y() - virtualMin.y());

    if(pi.x() > virtualMax.x() - fringe && pi.y() < virtualMin.y() + fringe)
      emitGhost(i, virtualMin.x() - virtualMax.x(), virtualMax.y() - virtualMin.y());

    if(pi.x() > virtualMax.x() - fringe && pi.y() > virtualMax.y() - fringe)
      emitGhost(i, virtualMin.x() - virtualMax.x(), virtualMin.y() - virtualMax.y());

    if(pi.x() < virtualMin.x() + fringe && pi.y() > virtualMax.y() - fringe)
      emitGhost(i, virtualMax.x() - virtualMin.x(), virtualMin.y() - virtualMax.y());

  }
}

// std::cout << "Removed " << filteredCount << " ptcls, added " << ghostCount << "ghost ptcls, diff " << filteredCount - ghostCount << " actual diff " << oldGhostCount - filteredCount << std::endl;
// oldGhostCount = ghostCount;
return;

// std::cout << "---- Beginning Ghost Generation ----\n";
auto mapCells = [&numPtcls, this](int32_t ci, int32_t cj, int32_t cti, int32_t ctj, bool flipX = false, bool flipY = false){
    vec cpos{ci * baseSupport, cj * baseSupport};
    vec tpos{cti * baseSupport, ctj * baseSupport};

    auto& cell = getCell(ci, cj);
    for(auto srcIdx: cell){
      vec relPos = fluidPosition[srcIdx] - cpos;
      vec relative_position = relPos / baseSupport;
      // if(flipX)
      //   relative_position.x() = 1. - relative_position.x();
      // if(flipY)
      //   relative_position.y() = 1. - relative_position.y();
      relPos = relative_position * baseSupport;

      fluidPosition[numPtcls] = relPos + tpos;
      // printf("Emitting Ghost %05d for %05d from [%g %g] to [%g %g]\n", numPtcls,srcIdx,fluidPosition[srcIdx].x(), fluidPosition[srcIdx].y(),fluidPosition[numPtcls].x(),fluidPosition[numPtcls].y());
      fluidVelocity[numPtcls] = fluidVelocity[srcIdx];
      fluidAngularVelocity[numPtcls] = fluidAngularVelocity[srcIdx];
      fluidVorticity[numPtcls] = fluidVorticity[srcIdx];
      fluidPriorPressure[numPtcls] = fluidPriorPressure[srcIdx];
      fluidArea[numPtcls] = fluidArea[srcIdx];
      fluidRestDensity[numPtcls] = fluidRestDensity[srcIdx];
      fluidUID[numPtcls]=fluidUID[srcIdx];
      fluidInitialPosition[numPtcls] = fluidPosition[numPtcls];
      fluidSupport[numPtcls] = fluidSupport[srcIdx];
      fluidGhostIndex[numPtcls] = srcIdx;
      ghostParticles.push_back(numPtcls);
      numPtcls++;
    }

};

// map cell B
if(periodicY){
  for(int32_t ci = 1; ci < cellsX-1; ++ci)
  for(int32_t cj = cellsY-2;cj < cellsY-1; ++cj)
  mapCells(ci,cj, ci, 0,false,true);
}
// map cell D
if(periodicX){
  for(int32_t ci = cellsX-2; ci < cellsX-1; ++ci)
  for(int32_t cj = 1;cj < cellsY-1; ++cj)
  mapCells(ci,cj, 0, cj,true,false);
}
// map cell F
if(periodicY){
  for(int32_t ci = 1; ci < cellsX-1; ++ci)
  for(int32_t cj = 1;cj < 2; ++cj)
  mapCells(ci,cj, ci, cellsY-1,false,true);
}
// map cell H
if(periodicX){
  for(int32_t ci = 1; ci < 2; ++ci)
  for(int32_t cj = 1;cj < cellsY-1; ++cj)
  mapCells(ci,cj, cellsX-1, cj,true,false);
}
if(periodicX && periodicY){
  mapCells(1,1,cellsX-1,cellsY-1,true);
  mapCells(cellsX-2,1,0,cellsY-1,true);
  mapCells(1,cellsY-2,cellsX-1,0,true);
  mapCells(cellsX-2,cellsY-2,0,0,true);
}
for(auto idx : ghostParticles){
  // std::cout << idx << " [ " << fluidGhostIndex[idx] << " ]\t";
}
// std::cout << std::endl;
for(int32_t c = 0; c < cellsX * cellsY; ++c){
    cellArray[c] = std::vector<int32_t>{};
}
}



fs::path expand(fs::path in) {
#ifndef _WIN32
  if (in.string().size() < 1)
    return in;

  const char *home = getenv("HOME");
  if (home == NULL) {
    std::cerr << "error: HOME variable not set." << std::endl;
    throw std::invalid_argument("error: HOME environment variable not set.");
  }

  std::string s = in.string();
  if (s[0] == '~') {
    s = std::string(home) + s.substr(1, s.size() - 1);
    return fs::path(s);
  } else {
    return in;
  }
#else
  if (in.string().size() < 1)
    return in;

  const char *home = getenv("USERPROFILE");
  if (home == NULL) {
    std::cerr << "error: USERPROFILE variable not set." << std::endl;
    throw std::invalid_argument("error: USERPROFILE environment variable not set.");
  }

  std::string s = in.string();
  if (s[0] == '~') {
    s = std::string(home) + s.substr(1, s.size() - 1);
    return fs::path(s);
  } else {
    return in;
  }
#endif
}

fs::path SPHSimulation::resolveFile(std::string fileName, std::vector<std::string> search_paths){
  fs::path working_dir = pm.get<std::string>("internal.working_directory");
  fs::path source_dir = pm.get<std::string>("internal.source_directory");
  fs::path build_dir = pm.get<std::string>("internal.build_directory");
  fs::path expanded = expand(fs::path(fileName));

  fs::path base_path = "";
  if (fs::exists(expand(fs::path(fileName))))
    return expand(fs::path(fileName));
  for (const auto &path : search_paths) {
    auto p = expand(fs::path(path));
    if (fs::exists(p / fileName))
      return p.string() + std::string("/") + fileName;
  }

  if (fs::exists(fileName))
    return fs::path(fileName);
  if (fs::exists(expanded))
    return expanded;

  for (const auto &pathi : search_paths) {
    auto path = expand(fs::path(pathi));
    if (fs::exists(working_dir / path / fileName))
      return (working_dir / path / fileName).string();
    if (fs::exists(source_dir / path / fileName))
      return (source_dir / path / fileName).string();
    if (fs::exists(build_dir / path / fileName))
      return (build_dir / path / fileName).string();
  }

  if (fs::exists(working_dir / fileName))
    return (working_dir / fileName);
  if (fs::exists(source_dir / fileName))
    return (source_dir / fileName);
  if (fs::exists(build_dir / fileName))
    return (build_dir / fileName);

  std::stringstream sstream;
  sstream << "File '" << fileName << "' could not be found in any provided search path" << std::endl;
  std::cerr << sstream.str();
  throw std::runtime_error(sstream.str().c_str());
}