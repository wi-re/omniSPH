#include <simulation/SPH.h>


void SPHSimulation::emitParticles(){
    auto currTime = pm.get<double>("sim.time");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& dt = pm.get<scalar>("sim.dt");
    currTime = currTime - dt;

        auto triangleSign = [](vec v1, vec v2, vec v3){
            return (v1.x() - v3.x()) * (v2.y() - v3.y()) - (v2.x() - v3.x()) * (v1.y() - v3.y());
        };
        auto pointInTriangle = [triangleSign](vec pt, vec v1, vec v2, vec v3){
            auto d1 = triangleSign(pt, v1, v2);
            auto d2 = triangleSign(pt, v2, v3);
            auto d3 = triangleSign(pt, v3, v1);

            bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

            return !(has_neg && has_pos);
        };

    for(auto& source: fluidSources){

        if(source.emitter != emitter_t::inlet) continue;
        if(source.timeLimit > 0. && currTime > source.timeLimit) continue;


        // auto speed = source.emitterVelocity.norm();
        auto velDir = source.emitterVelocity.normalized();
        auto speed = source.emitterRampTime > 0. ? std::clamp(currTime / source.emitterRampTime, 0., 1.) * source.emitterVelocity.norm() : source.emitterVelocity.norm();

        if(inletOnce)
            source.emitterOffset = -dt * speed;
        auto support = std::sqrt(source.emitterRadius * source.emitterRadius  * targetNeighbors);
        source.emitterOffset = std::fmod(source.emitterOffset + dt * speed, 1.0 * packing_2D * support);
        // source.emitterOffset = source.emitterOffset + dt * speed;
        // printf("Offset: %g, Period: %g, Increment: %g\n", source.emitterOffset, 4. * packing_2D * support, dt *speed);

        auto ptcls = source.genParticles();
        for(int32_t i = 0; i < ptcls.size(); ++ i){
            auto pos = ptcls[i];
            pos += velDir * source.emitterOffset;
        auto [ix, iy] = getCellIdx(ptcls[i].x(), ptcls[i].y());
        bool emit = true;
        for (int32_t xi = -2; xi <= 2; ++xi) {
            for (int32_t yi = -2; yi <= 2; ++yi) {
                if(ix +xi < 0 || ix +xi >= cellsX || iy+yi <0 || iy + yi >= cellsY) continue;
                const auto& cell = getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    auto& pj = fluidPosition[j];
                    vec r = pj - pos;
                    if (r.squaredNorm() <= /*1.5 * 2.0 * 2.0 **/ 0.5* packing_2D * packing_2D * support * support) {
                        emit = false; 
                    }
                }
            }
        }
        if(pos.x() <= source.emitterMin.x() || pos.x() >= source.emitterMax.x() || pos.y() <= source.emitterMin.y() || pos.y() >= source.emitterMax.y())
        emit = false;

        for(auto t: boundaryTriangles)
            if(pointInTriangle(pos, t.v0, t.v1, t.v2))
            emit = false;
        

        if(!emit) continue;
            fluidPosition[numPtcls]         = pos;
            fluidVelocity[numPtcls]         = source.emitterVelocity;
            fluidArea[numPtcls]             = source.emitterRadius * source.emitterRadius * double_pi;
            fluidRestDensity[numPtcls]      = source.emitterDensity;
            fluidSupport[numPtcls]          = std::sqrt(fluidArea[numPtcls] * targetNeighbors / double_pi);
            fluidPriorPressure[numPtcls]    = 0.;
            fluidVorticity[numPtcls] = 0.;
            fluidAngularVelocity[numPtcls] = 0.;
            fluidUID[numPtcls] = fluidCounter++;
            fluidInitialPosition[numPtcls] = pos;
            fluidGhostIndex[numPtcls] = -1;

            // getCell(pos.x(), pos.y()).push_back(numPtcls);

            numPtcls += 1;
        }
    }



    inletOnce = false;
}


std::vector<vec> fluidSource::genParticles() const{
    scalar area = double_pi * emitterRadius * emitterRadius;
    scalar support = std::sqrt(area * targetNeighbors / double_pi);
    scalar packing = packing_2D * support / compressionRatio;
    printf("Compression: %g | %g | %g\n", compressionRatio, packing, packing_2D * support);
    printf("radius: %g, area: %g, support: %g, packing: %g\n", emitterRadius, area, support, packing);
    vec center = (emitterMin + emitterMax) / 2.;
    vec pdiff = (emitterMax - emitterMin) / 2.;
    scalar radius = std::min(pdiff.x(), pdiff.y());

  auto gen_position = [](auto r, auto i, auto j) -> vec {
    //return vec(r * (2.0 * scalar(i) + scalar(j % 2)), r * ::sqrt(3.0) * scalar(j));
    return vec(r * scalar(i) , r * scalar(j));
  };
  auto diff = emitterMax - emitterMin;
  auto requiredSlices_x = ::ceil(diff.x() / packing);
  //auto requiredSlices_y = ::ceil(diff.y() / (::sqrt(3.0) * packing));
  auto requiredSlices_y = ::ceil(diff.y() / packing);
 std::cout << "Generating particles on a " << requiredSlices_x << " x " << requiredSlices_y << " grid " << std::endl;
  std::vector<vec> points;
  for (int32_t x_it = 0; x_it < requiredSlices_x + 1; ++x_it)
    for (int32_t y_it = 0; y_it < requiredSlices_y + 1; ++y_it) {
      vec p = emitterMin;
      vec g = gen_position(packing, x_it, y_it);
      vec pos = p + g;
      if(shape == shape_t::spherical){
                auto d = pos - center;
                auto l = d.norm();
                //std::cout << pos.x() << " " << pos.y() << " -> " << d.x() << " " << d.y() << " -> " << l << " / " << radius << " -> " <<(l>= radius) << std::endl;
                if (l >= radius) continue;
                }
      if (pos.x() <= emitterMax.x() && pos.y() <= emitterMax.y())
        points.emplace_back(pos.x(), pos.y());
    }
  return points;
}

std::vector<vec> SPHSimulation::genParticles(vec minCoord, vec maxCoord, scalar radius, scalar packing){
    scalar area = double_pi * radius * radius;
    scalar support = std::sqrt(area * targetNeighbors / double_pi);

  auto gen_position = [](auto r, auto i, auto j) -> vec {
    //return vec(r * (2.0 * scalar(i) + scalar(j % 2)), r * ::sqrt(3.0) * scalar(j));
    return vec(r * scalar(i) , r * scalar(j));
  };
  packing *= support;
  //packing *= 0.95;
  //packing = support;
  auto diff = maxCoord - minCoord;
  auto requiredSlices_x = ::ceil(diff.x() / packing);
  //auto requiredSlices_y = ::ceil(diff.y() / (::sqrt(3.0) * packing));
  auto requiredSlices_y = ::ceil(diff.y() / packing);
 std::cout << "Generating particles on a " << requiredSlices_x << " x " << requiredSlices_y << " grid " << std::endl;
  std::vector<vec> points;
  for (int32_t x_it = 0; x_it < requiredSlices_x + 1; ++x_it)
    for (int32_t y_it = 0; y_it < requiredSlices_y + 1; ++y_it) {
      vec p = minCoord;
      vec g = gen_position(packing, x_it, y_it);
      vec pos = p + g;
      if (pos.x() < maxCoord.x() && pos.y() <= maxCoord.y()+ support * 0.2)
        points.emplace_back(pos.x(), pos.y());
    }
  return points;
}
std::vector<vec> SPHSimulation::genParticlesCircle(vec minCoord, vec maxCoord, scalar rad, scalar packing) {
    scalar area = double_pi * rad * rad;
    scalar support = std::sqrt(area * targetNeighbors / double_pi);
    vec center = (minCoord + maxCoord) / 2.;
    vec pdiff = (maxCoord - minCoord) / 2.;
    scalar radius = std::min(pdiff.x(), pdiff.y());

    auto gen_position = [](auto r, auto i, auto j) -> vec {
        //return vec(r * (2.0 * scalar(i) + scalar(j % 2)), r * ::sqrt(3.0) * scalar(j));
        return vec(r * scalar(i), r * scalar(j));
    };
    packing *= support;
    //packing *= 0.95;
    //packing = support;
    auto diff = maxCoord - minCoord;
    auto requiredSlices_x = ::ceil(diff.x() / packing);
    //auto requiredSlices_y = ::ceil(diff.y() / (::sqrt(3.0) * packing));
    auto requiredSlices_y = ::ceil(diff.y() / packing);
    // std::cout << "Generating particles on a " << requiredSlices_x << " x " << requiredSlices_y << " hex grid " << std::endl;
    std::vector<vec> points;
    for (int32_t x_it = 0; x_it < requiredSlices_x + 1; ++x_it)
        for (int32_t y_it = 0; y_it < requiredSlices_y + 1; ++y_it) {
            vec p = minCoord;
            vec g = gen_position(packing, x_it, y_it);
            vec pos = p + g;
            if (pos.x() < maxCoord.x() && pos.y() <= maxCoord.y() + support * 0.2) {
                auto d = pos - center;
                auto l = d.norm();
                //std::cout << pos.x() << " " << pos.y() << " -> " << d.x() << " " << d.y() << " -> " << l << " / " << radius << " -> " <<(l>= radius) << std::endl;
                if (l >= radius) continue;
                points.emplace_back(pos.x(), pos.y());
            }
        }
    return points;
}

std::vector<int32_t>& SPHSimulation::getCell(scalar x, scalar y){
    if (x != x)
        x = (scalar)0.0;
    if (y != y)
        y = (scalar)0.0;
    std::size_t xi = static_cast<std::size_t>(::floor(std::clamp(x - domainMin.x(), (scalar)0.0, domainMax.x() - domainMin.x()) / baseSupport));
    std::size_t yi = static_cast<std::size_t>(::floor(std::clamp(y - domainMin.y(), (scalar)0.0, domainMax.x() - domainMin.y()) / baseSupport));
    xi = std::clamp(xi, (std::size_t)0, cellsX - 1);
    yi = std::clamp(yi, (std::size_t)0, cellsY - 1);
    return cellArray[yi * (cellsX)+xi];
    }
std::pair<int32_t, int32_t> SPHSimulation::getCellIdx(scalar x, scalar y){
    if (x != x)
        x = (scalar)0.0;
    if (y != y)
        y = (scalar)0.0;
    std::size_t xi = static_cast<std::size_t>(::floor(std::clamp(x - domainMin.x(), (scalar)0.0, domainMax.x() - domainMin.x()) / baseSupport));
    std::size_t yi = static_cast<std::size_t>(::floor(std::clamp(y - domainMin.y(), (scalar)0.0, domainMax.x() - domainMin.y()) / baseSupport));
    xi = std::clamp(xi, (std::size_t)0, cellsX - 1);
    yi = std::clamp(yi, (std::size_t)0, cellsY - 1);
    return std::make_pair((int32_t)xi, (int32_t)yi);}
std::vector<int32_t>& SPHSimulation::getCell(int32_t xi, int32_t yi){
    return cellArray[yi * (cellsX)+xi];
    }
std::vector<int32_t>& SPHSimulation::getBoundaryCell(int32_t xi, int32_t yi){
    return cellBoundaryArray[yi * (cellsX)+xi];
    }
std::vector<int32_t>& SPHSimulation::getTriangleCell(int32_t xi, int32_t yi){
    return cellTriangleArray[yi * (cellsX)+xi];
    }
