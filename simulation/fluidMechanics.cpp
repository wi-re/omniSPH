#include <simulation/SPH.h>
#include <simulation/2DMath.h>
#include <cfloat>

void SPHSimulation::density(){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        auto& pi = fluidPosition[i];
        auto& hi = fluidSupport[i];
        auto rho = 0.;
        for (int32_t j : fluidNeighborList[i]) {
            rho += fluidArea[j] * W(pi, fluidPosition[j], hi, fluidSupport[j]);
        }
        //pi.rho = std::max(pi.rho, 0.5);
        boundaryFunc(i, [&rho](auto bpos, auto d, auto k, auto gk, auto triangle) {
            rho += k; });

        fluidDensity[i] = rho;
    }
}
void SPHSimulation::Integrate(){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& dt = pm.get<scalar>("sim.dt");


        auto& time = pm.get<scalar>("sim.time");
        auto delay = 0.13;
        //if (time < delay)
        //    speed2 = 0.0;
        //else
        //else if (time < 1.0 + delay) {
        //    speed2 = speedGoal * (time - delay);
        //}
    
        //speed2 = speed;
        auto spacing_2D = 0.23999418487168855;

        for(auto& source : fluidSources){
            if(source.emitter != emitter_t::velocitySource) continue;
            
        if(source.timeLimit > 0. && time > source.timeLimit) continue;
    vec center = (source.emitterMin + source.emitterMax) / 2.;
    vec pdiff = (source.emitterMax - source.emitterMin) / 2.;
    scalar radius = std::min(pdiff.x(), pdiff.y());

            for (int32_t i = 0; i < numPtcls; ++i) {
                auto pos = fluidPosition[i];
                auto speed = source.emitterRampTime > 0. ? std::clamp(time / source.emitterRampTime, 0., 1.) * source.emitterVelocity.norm() : source.emitterVelocity.norm();

                if(source.shape == shape_t::rectangular)
                if(pos.x() <= source.emitterMax.x() && pos.x() >= source.emitterMin.x() && pos.y() <= source.emitterMax.y() && pos.y() >= source.emitterMin.y()){
                    auto mu = 3.5;
                    auto xr = (pos.x() - source.emitterMin.x()) / (source.emitterMax.x() - source.emitterMin.x());
                    if(source.emitterMin.x() > .5)
                    xr = 1. - xr;
                    auto gamma = 1. - (std::exp(std::pow(xr,mu))-1.) / (std::exp(1) - 1.);

                    auto emitterVel = source.emitterVelocity;
                    fluidVelocity[i] += dt * fluidAccel[i];
                    fluidVelocity[i] = fluidVelocity[i] * ( 1. - gamma) + gamma * emitterVel;
                    // printf("%g -> %g : %g -> %g\n", pos.x(), xr, std::exp(std::pow(xr,mu)), gamma );
                    fluidAngularVelocity[i] = 0.;
                    fluidAccel[i] = vec(0.,0.);
                }
                if(source.shape == shape_t::spherical){
                auto d = pos - center;
                auto l = d.norm();
                //std::cout << pos.x() << " " << pos.y() << " -> " << d.x() << " " << d.y() << " -> " << l << " / " << radius << " -> " <<(l>= radius) << std::endl;
                if (l >= radius) continue;
                    fluidVelocity[i] = source.emitterVelocity.normalized() * speed;
                    fluidAngularVelocity[i] = 0.;
                    fluidAccel[i] = vec(0.,0.);
                }
            }

        }

        for(auto& source : fluidSources){
            if(source.emitter != emitter_t::inlet) continue;
            
        if(source.timeLimit > 0. && time > source.timeLimit) continue;
    vec center = (source.emitterMin + source.emitterMax) / 2.;
    vec pdiff = (source.emitterMax - source.emitterMin) / 2.;
    scalar radius = std::min(pdiff.x(), pdiff.y());

            for (int32_t i = 0; i < numPtcls; ++i) {
                auto pos = fluidPosition[i];
                auto speed = source.emitterRampTime > 0. ? std::clamp(time / source.emitterRampTime, 0., 1.) * source.emitterVelocity.norm() : source.emitterVelocity.norm();

                if(source.shape == shape_t::rectangular)
                if(pos.x() <= source.emitterMax.x() && pos.x() >= source.emitterMin.x() && pos.y() <= source.emitterMax.y() && pos.y() >= source.emitterMin.y()){
                    auto emitterVel = source.emitterVelocity;
                    fluidVelocity[i] = emitterVel;
                    // printf("%g -> %g : %g -> %g\n", pos.x(), xr, std::exp(std::pow(xr,mu)), gamma );
                    fluidAngularVelocity[i] = 0.;
                    fluidAccel[i] = vec(0.,0.);
                }
                if(source.shape == shape_t::spherical){
                auto d = pos - center;
                auto l = d.norm();
                //std::cout << pos.x() << " " << pos.y() << " -> " << d.x() << " " << d.y() << " -> " << l << " / " << radius << " -> " <<(l>= radius) << std::endl;
                if (l >= radius) continue;
                    fluidVelocity[i] = source.emitterVelocity.normalized() * speed;
                    fluidAngularVelocity[i] = 0.;
                    fluidAccel[i] = vec(0.,0.);
                }
            }

        }



        auto& damping = pm.get<scalar>("props.damping");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        auto& pi = fluidPosition[i];
        auto& vi = fluidVelocity[i];
        auto& ai = fluidAccel[i];
        vi += dt * ai;
        vi *= (1.0 - damping);
        pi += dt * vi;
    }

        scalar vMax = 0.0;
        scalar minH = DBL_MAX;
        //auto& dt = pm.get<scalar>("sim.dt");
        for (int32_t i = 0; i < numPtcls; ++i) {
            vMax = std::max(vMax,fluidVelocity[i].norm());
            minH = std::min(minH, fluidSupport[i]);
        }
    

        auto& t = pm.get<scalar>("sim.time");
        t += dt;
        pm.get<int32_t>("sim.frame")++;
        auto& vMaxr = pm.get<scalar>("ptcl.maxVelocity");
        vMaxr = vMax;
        auto& dtmin = pm.get<scalar>("sim.minDt");
        auto& dtmax = pm.get<scalar>("sim.maxDt");
        dt = std::clamp(0.4 * minH / vMax, dtmin, dtmax);
}


#define rhoj fluidRestDensity[i]
// #define rhoj fluidRestDensity[j]

void SPHSimulation::XSPH(){
     auto& numPtcls = pm.get<int32_t>("props.numPtcls");
     auto& viscosityConstant = pm.get<scalar>("ptcl.viscosityConstant");
     auto& boundaryViscosity = pm.get<scalar>("ptcl.boundaryViscosity");
     std::vector<vec> tempV;
     //#pragma omp parallel for
     for (int32_t i = 0; i < numPtcls; ++i) {
         auto& pi = fluidPosition[i];
         auto& vi = fluidVelocity[i];
         auto& hi = fluidSupport[i];
         tempV.push_back(fluidVelocity[i]);
         for (auto& j : fluidNeighborList[i]) {
             auto& pj = fluidPosition[j];
             auto& hj = fluidSupport[j];
             tempV[i] += viscosityConstant * (fluidDensity[i] + fluidDensity[j]) * fluidArea[j] / 
                 (fluidDensity[i] + fluidDensity[j]) * scalar(2.0) * W(pi, pj, hi, hj) * (fluidVelocity[j] - fluidVelocity[i]);
         }
     }

 #pragma omp parallel for
     for (int32_t i = 0; i < numPtcls; ++i)
         fluidVelocity[i] = tempV[i];


}
void SPHSimulation::computeAlpha(bool density){
    auto& dt = pm.get<scalar>("sim.dt");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        vec kernelSum1(0, 0);
        scalar kernelSum2 = 0.0;
        if (density)
            boundaryFunc(i, [&kernelSum1, &kernelSum2](auto bpos, auto d, auto k, auto gk, auto triangle) { kernelSum1 += gk; });
        for (int32_t j : fluidNeighborList[i]) {
            vec kernel = gradW(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]);
            kernelSum1 += fluidActualArea[j] * kernel;
            kernelSum2 += fluidActualArea[j] * fluidActualArea[j] / (fluidArea[j] * rhoj) * kernel.dot(kernel);
        }
        fluidAlpha[i] = -dt * dt * fluidActualArea[i] / (fluidArea[i] * fluidRestDensity[i]) * kernelSum1.dot(kernelSum1) - dt * dt * fluidActualArea[i] * kernelSum2;
   
    }
}
void SPHSimulation::computeSourceTerm(bool density){
    auto& dt = pm.get<scalar>("sim.dt");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        scalar sourceTerm = density ? scalar(1.0) - fluidDensity[i] : scalar(0.0);
        if (density)
            boundaryFunc(i, [this,i, dt, &sourceTerm](auto bpos, auto d, auto k, auto gk, auto triangle) {
            sourceTerm = sourceTerm - dt * fluidPredVelocity[i].dot(gk);
                });
        for (int32_t j : fluidNeighborList[i]) {
            sourceTerm -= dt * fluidActualArea[j] * (fluidPredVelocity[i] - fluidPredVelocity[j]).dot(gradW(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]));
        }
        fluidSourceTerm[i] = sourceTerm;
        fluidPressure2[i] = (scalar)0.0;
    }
}
void SPHSimulation::computeAcceleration(bool density){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        vec kernelSum((scalar)0.0, (scalar)0.0);
        for (int32_t j : fluidNeighborList[i]) {
            kernelSum += -fluidArea[j]* rhoj * (fluidPressure2[i] / power(fluidDensity[i] * fluidRestDensity[i], 2) + fluidPressure2[j] / power(fluidDensity[j] * rhoj, 2)) *
                gradW(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]);
        }
        fluidPredAccel[i] += kernelSum;
        fluidPressure1[i] = fluidPressure2[i];
    }
}
void SPHSimulation::updatePressure(bool density ){
    auto& dt = pm.get<scalar>("sim.dt");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& backgroundPressureSwitch = pm.get<bool>("props.backgroundPressure");
    auto& vMaxr = pm.get<scalar>("ptcl.maxVelocity");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        scalar kernelSum = (scalar)0.0;
        if (density)
            boundaryFunc(i,
                [this,dt, i, &kernelSum](auto bpos, auto d, auto k, auto gk, auto triangle) { kernelSum += dt * dt * fluidPredAccel[i].dot(gk); });
        for (int32_t j : fluidNeighborList[i]) {
            kernelSum += dt * dt * fluidActualArea[j] * (fluidPredAccel[i] - fluidPredAccel[j]).dot(gradW(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]));
        }
        scalar omega = (scalar)0.5;
        scalar pressure = fluidPressure1[i] + omega / fluidAlpha[i] * (fluidSourceTerm[i] - kernelSum);
        //if(di.pressureBoundary == 0.0)
        if (backgroundPressureSwitch && density)
            pressure += 0.5 * vMaxr * vMaxr * fluidRestDensity[i];
        pressure = density ? std::max(pressure, (scalar)0.0) : pressure;
        scalar residual = kernelSum - fluidSourceTerm[i];
        if (::abs(fluidAlpha[i]) < 1e-25 || pressure != pressure || pressure > 1e25)
            pressure = residual = 0.0;
        //pressure = std::max(pressure, speed * speed * rho0);
        fluidPressure2[i] = pressure;
        fluidDpDt[i] = std::max(residual, (scalar)-0.001) * fluidArea[i];
        fluidDensityStar[i] = residual * fluidArea[i];
    }
}
scalar SPHSimulation::calculateBoundaryPressureMLS(int32_t i, vec pb, bool density) {
    vec vecSum(0, 0), d_bar(0, 0);
    matrix M = matrix::Zero();
    scalar sumA = (scalar)0.0, sumB = (scalar)0.0, d_sum = (scalar)0.0;
    auto support = fluidSupport[i];
    int32_t ii = 0;
    auto [ix, iy] = getCellIdx(pb.x(), pb.y());
    for (int32_t xi = -1; xi <= 1; ++xi) {
        for (int32_t yi = -1; yi <= 1; ++yi) {
            if (xi + ix < 0 || xi + ix >= cellsX) continue;
            if (yi + iy < 0 || yi + iy >= cellsY) continue;
            const auto& cell = getCell(ix + xi, iy + yi);
            for (auto j : cell) {
                vec r = fluidPosition[j] - pb;
                if (r.squaredNorm() <= support * support) {
                    ii++;
                    scalar fac = fluidArea[j] * W(fluidPosition[j], pb, support, support);
                    d_bar += fluidPosition[j] * fac;
                    d_sum += fac;
                }
            }
        }
    }if (ii == 0) return 0.;
    d_bar /= d_sum;
    vec x_b = pb - d_bar;
    //std::cout << "---------------------\n";
    for (int32_t xi = -1; xi <= 1; ++xi) {
        for (int32_t yi = -1; yi <= 1; ++yi) {
            if (xi + ix < 0 || xi + ix >= cellsX) continue;
            if (yi + iy < 0 || yi + iy >= cellsY) continue;
            const auto& cell = getCell(ix + xi, iy + yi);
            for (auto j : cell) {
                vec r = fluidPosition[j] - pb;
                if (r.squaredNorm() <= (scalar)support * support) {
                    scalar Wbbf = W(fluidPosition[j], pb,support,support);
                    vec pjb = fluidPosition[j] - d_bar;
                    M += (pjb) * (pjb).transpose();
                    vecSum += pjb * fluidPressure2[j] * fluidArea[j] * Wbbf;
                    sumA += fluidPressure2[j] * fluidArea[j] * Wbbf;
                    sumB += fluidArea[j] * Wbbf;
                }
            }
        }
    }
    auto [U, S, V] = svd2x2(M);
    S(0, 0) = abs(S(0, 0)) > (scalar)1e-2f ? (scalar)1.0 / S(0, 0) : (scalar)0.0;
    S(1, 1) = abs(S(1, 1)) > (scalar)1e-2f ? (scalar)1.0 / S(1, 1) : (scalar)0.0;
    S(0, 1) = S(1, 0) = (scalar)0.0;
    matrix Mp = V * S * U.transpose();
    scalar alpha = sumA / sumB;
    auto beta = Mp(0, 0) * vecSum.x() + Mp(0, 1) * vecSum.y();
    auto gamma = Mp(1, 0) * vecSum.x() + Mp(1, 1) * vecSum.y();
    scalar det = M.determinant();
    //std::cout << std::defaultfloat;
    // std::cout << alpha << " " << beta << " " << gamma << " [ " << Mp(0, 0) << " " << Mp(0, 1) << " ] [ " << Mp(1, 0) << " " << Mp(1, 1) << " ] " << " " << det << " " << x_b.x() << " " << x_b.y() << "\n";
    if (det != det)
        beta = gamma = (scalar)0.0;
    auto pressure = alpha + beta * x_b.x() + gamma * x_b.y();
    if (pressure != pressure || alpha != alpha)
        pressure = (scalar)0.0;
    pressure = density ? std::max(pressure, (scalar)0.0) : pressure;
    //if (pressure > 1e10)
      //  printParticle(i);
    return pressure;
}

void SPHSimulation::computeBoundaryTrianglePressure(bool density) {
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& baryPressure = pm.get<bool>("sim.barycentricPressure");
    if (!baryPressure) return;
#pragma omp parallel for
    for (int32_t i = 0; i < boundaryTriangles.size(); ++i) {
        const auto [t0, t1, t2] = boundaryTriangles[i];
        auto f0 = calculateBoundaryPressureMLS(-1, t0, true);
        auto f1 = calculateBoundaryPressureMLS(-1, t1, true);
        auto f2 = calculateBoundaryPressureMLS(-1, t2, true);
        boundaryBarycentricPressure[i] = std::make_tuple(f0, f1, f2);
    }

}
void SPHSimulation::computeBoundaryPressure(bool density ){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& baryPressure = pm.get<bool>("sim.barycentricPressure");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        fluidBoundaryPressure[i] = 0.0;
        fluidPredAccel[i] = vec(0, 0);
        if (!density)
            continue;

        if (baryPressure) {
            for (auto ti : fluidTriangleNeighborList[i]) {
                const auto& tri = boundaryTriangles[ti];
                auto [f0, f1, f2] = boundaryBarycentricPressure[ti];
                scalar fluidPressure = density ? std::max((scalar)0.0, fluidPressure2[i]) : fluidPressure2[i];
                auto [pb, d] = closestPointTriangle(fluidPosition[i], tri);
                scalar pressure = calculateBoundaryPressureMLS(0, pb, density);
                //f0 = f1 = f2 = pressure;
                //f0 = std::max(f0,std::max(f1, f2));
                //f1 = f2 = f0;
                auto [hit, pb2, d2, k, gk] = interactTriangleBaryCentric(fluidPosition[i], fluidSupport[i], fluidRestDensity[i], tri, fluidDensity[i] * fluidRestDensity[i], fluidPressure, f0, f1, f2);
                //auto [hit2, pb3, d3, k3, gk3] = interactTriangle(fluidPosition[i], tri);


                if (hit) {
                    //auto [hit, pb, d, k2, gk2] = interactTriangleBaryCentric(fluidPosition[i], tri, fluidDensity[i] * rho0, fluidPressure, pressure, pressure, pressure);
                    //std::cout << "Hit Triangle!!!!\n\n\n\n";
                    //if (std::abs(gk.norm() - gk2.norm())>1e-2)
                    //{
                    //    std::cout << "\n";
                    //    std::cout << std::setprecision(5) << f0 << " " << f1 << " " << f2 << " -> " << pressure << " => " << k << " @ [ " << gk.x() << " " << gk.y() << "] : " << k3 << " @ [ " << gk3.x() << " " << gk3.y() << "] : ";
                    //    std::cout << k2 << " @ [ " << gk2.x() << " " << gk2.y() << "] \n";
                    //}
                    fluidPredAccel[i] += -scalar(1.0) * gk;
                }
                //if (f0 > 1e12 || f1 > 1e12 || f2 > 1e12) {
                //    std::cout << std::defaultfloat;
                //    std::cout << "\n\n";
                //    std::cout << f0 << " " << f1 << " " << f2 << " - " << fluidPressure << std::endl;
                //    std::cout << fluidPosition[i].x() << " : " << fluidPosition[i].y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << d << ", " << k << " ==> " << gk.x() << " : " << gk.y() << std::endl;
                //    printParticle(i);
                //}

                fluidBoundaryPressure[i] += f0 + f1 + f2;
            }
        }
        else {
            boundaryFunc(i, [&, this](auto pb, auto d, auto k, auto gk, auto triangle) {
                scalar pressure = calculateBoundaryPressureMLS(i, pb, density);
                scalar fluidPressure = density ? std::max((scalar)0.0, fluidPressure2[i]) : fluidPressure2[i];
                scalar boundaryPressure = density ? std::max((scalar)0.0, pressure) : pressure;
                // boundaryPressure = fluidPressure;
                fluidBoundaryPressure[i] += boundaryPressure;
                fluidPredAccel[i] +=
                    -scalar(1.0) * fluidRestDensity[i] * (fluidPressure / power(fluidDensity[i] * fluidRestDensity[i], 2) + boundaryPressure / (fluidRestDensity[i] * fluidRestDensity[i])) * gk;
                });
        }
    }
}
void SPHSimulation::predictVelocity(bool density ){
    auto& dt = pm.get<scalar>("sim.dt");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        fluidPredVelocity[i] = fluidVelocity[i] + dt * fluidAccel[i];
        fluidActualArea[i] = fluidArea[i] / fluidDensity[i];
    }
}
void SPHSimulation::updateVelocity(bool density ){
    auto& dt = pm.get<scalar>("sim.dt");
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
       fluidAccel[i] += fluidPredAccel[i];
        fluidPredVelocity[i] += dt * fluidPredAccel[i];
    }
}
int32_t SPHSimulation::divergenceSolve() {
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& active = pm.get<bool>("dfsph.divergenceSolve");
    if (!active) return 0;
    //return 0;
    predictVelocity();
    computeAlpha(false);
    computeSourceTerm(false);
    //scalar error = (scalar)0.0;
    //int32_t counter = 0;
    auto& limit = pm.get<scalar>("dfsph.divergenceEta");
    auto& error = pm.get<scalar>("dfsph.divergenceError");
    auto& counter = pm.get<int32_t>("dfsph.divergenceIterations");
    counter = 0;
    do {
        //computeBoundaryTrianglePressure(false);
        computeBoundaryPressure(false);
        computeAcceleration(false);
        updatePressure(false);
        error = (scalar)0.0;
        for (int32_t i = 0; i < numPtcls; ++i)
            error += fluidDpDt[i] / fluidArea[i];
        error /= numPtcls;
        //std::cout << "Divergence: " << counter << " -> " << error << std::endl;
    } while (counter++ < 3 || (error > (scalar)limit && counter < 3));
    computeBoundaryPressure(false);
    computeAcceleration(false);
    updateVelocity(false);
    return counter;
}
int32_t SPHSimulation::densitySolve(){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    // return PCISPH();
    predictVelocity();
    computeAlpha(true);
    computeSourceTerm(true);
    //scalar error = (scalar)0.0;
    //int32_t counter = 0;
    for (int32_t i = 0; i < numPtcls; ++i){
       fluidPressure1[i] = fluidPressure2[i] = 0.5 * fluidPressure1[i];
    }
    auto& limit = pm.get<scalar>("dfsph.densityEta");
    auto& error = pm.get<scalar>("dfsph.densityError");
    auto& counter = pm.get<int32_t>("dfsph.densityIterations");
    counter = 0;
    do {
        computeBoundaryTrianglePressure(true);
        computeBoundaryPressure(true);
        computeAcceleration(true);
        updatePressure(true);
        error = (scalar)0.0;
        for (int32_t i = 0; i < numPtcls; ++i)
            error += std::max(-0.001 * fluidArea[i], fluidDensityStar[i]) / fluidArea[i];
        error /= numPtcls;
        // std::cout << "Density: " << counter << " -> " << error << std::endl;
    } while (counter++ < 3 || (error > (scalar)limit && counter < 256));
    computeBoundaryPressure(true);
    computeAcceleration(true);
    updateVelocity();


    for (int32_t i = 0; i < numPtcls; ++i) {
        fluidPriorPressure[i] = fluidPressure1[i];
    }

    return counter;; }


void SPHSimulation::computeVorticity() {
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        scalar voriticity = 0.0;

        for (int32_t j : fluidNeighborList[i]) {
            auto grad = gradW(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]);
            auto vel = fluidVelocity[i] - fluidVelocity[j];
            auto term = fluidArea[j] / fluidDensity[j] * (vel.x() * grad.y() - vel.y() * grad.x());

            voriticity += term;
        }
        fluidVorticity[i] = voriticity;
    }
}
void SPHSimulation::refineVorticity() {
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    auto& nu_t = pm.get<scalar>("vorticity.nu_t");
    auto& intertiaInverse = pm.get<scalar>("vorticity.inverseInertia");
    auto& angularViscosity = pm.get<scalar>("vorticity.angularViscosity");

    std::vector<scalar> dwdt(numPtcls, 0.0);
    std::vector<vec> dvdt(numPtcls, vec(0.0, 0.0));
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        scalar voriticity = 0.0;
        vec vterm(0.0, 0.0);

        for (int32_t j :fluidNeighborList[i]) {
            auto grad = gradW(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]);
            auto vel = fluidVelocity[i] - fluidVelocity[j];
            auto vor = fluidAngularVelocity[i] - fluidAngularVelocity[j];
            auto term = fluidArea[j] / fluidDensity[j] * (vel.x() * grad.y() - vel.y() * grad.x());

            voriticity += term;
            vterm.x() += fluidArea[j] / fluidDensity[j] * (-vor * grad.y());
            vterm.y() += fluidArea[j] / fluidDensity[j] * (vor * grad.x());


        }
        dwdt[i] = nu_t * (voriticity - 2.0 * fluidAngularVelocity[i]) * intertiaInverse;
        dvdt[i] = nu_t * vterm;
    }

    auto& dt = pm.get<scalar>("sim.dt");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        fluidAngularVelocity[i] += dwdt[i] * dt;
        fluidVelocity[i] += dvdt[i] * dt;
    }

#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        scalar angularTerm = 0.0;
        for (int32_t j :fluidNeighborList[i]) {
            auto kernel = W(fluidPosition[i], fluidPosition[j], fluidSupport[i], fluidSupport[j]);
            auto vor = fluidAngularVelocity[j] - fluidAngularVelocity[i];
            auto term = fluidArea[j] / fluidDensity[j] * vor * kernel;
            angularTerm += term;
        }
        dwdt[i] = angularViscosity * angularTerm;
    }
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        fluidAngularVelocity[i] += dwdt[i];
    }
}
void SPHSimulation::externalForces() {
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");

    for(auto& gravity: gravitySources){
        if(gravity.pointSource){
            auto location = gravity.location;
            for (int32_t i = 0; i < numPtcls; ++i){
                auto dir = fluidPosition[i] - location;
                vec n = dir.normalized();
                fluidAccel[i] += -gravity.magnitude * n;
            }
        }
        else{
            for (int32_t i = 0; i < numPtcls; ++i)
                fluidAccel[i] += gravity.direction * gravity.magnitude;
        }
    }

}
void SPHSimulation::BXSPH() {
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");

    auto& viscosityConstant = pm.get<scalar>("ptcl.viscosityConstant");
    auto& boundaryViscosity = pm.get<scalar>("ptcl.boundaryViscosity");
    std::vector<vec> tempV;

    auto mu = boundaryViscosity;
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        if (fluidTriangleNeighborList[i].size() == 0) continue;
        auto f_visco = vec(0, 0);
        auto grad = vec(0, 0);
        auto ksum = 0.;
        auto vsum = vec(0, 0);
        for (auto ti : fluidTriangleNeighborList[i]) {
            const auto& t = boundaryTriangles[ti];
            auto [hit, pb2, d2, k, gk] = interactTriangle(fluidPosition[i], fluidSupport[i], t);
            if (hit) {
                grad += gk;
                ksum += k;
                vsum += gk;
            }
        }
        grad = -grad.normalized();

        auto v = fluidVelocity[i];
        auto proj = v.dot(grad) * grad;
        auto orth = v - proj;
        f_visco -= std::min(mu * ksum, 1.0) * orth;

        //f_visco += mu * orth;


        //for (auto ti : pi.neighborTriangles) {
        //    const auto& t = triangles[ti];

        //    const auto f0 = (pi.vel - vec(0, 0)).dot(pi.pos - t.v0) / ((pi.pos - t.v0).dot(pi.pos - t.v0) + 0.01 * support * support);
        //    const auto f1 = (pi.vel - vec(0, 0)).dot(pi.pos - t.v1) / ((pi.pos - t.v1).dot(pi.pos - t.v1) + 0.01 * support * support);
        //    const auto f2 = (pi.vel - vec(0, 0)).dot(pi.pos - t.v2) / ((pi.pos - t.v2).dot(pi.pos - t.v2) + 0.01 * support * support);

        //    auto [hit, pb2, d2, k, gk] = interactTriangle(pi.pos, t);

        //   // auto [hit3, pb32, d32, k2, gk2] = interactTriangle(pi.pos, t);
        //    if (hit) {
        //        //printf("%d: [%g %g] @ [%g %g] -> %g %g %g => [%g %g]\n", i, pi.pos.x(), pi.pos.y(), pi.vel.x(), pi.vel.y(), f0, f1, f2, gk.x(), gk.y());
        //        auto vc = (t.v0 + t.v1 + t.v2) / 3.;
        //        auto diff = pi.pos - vc;
        //        //auto grad = gk2.normalized();
        //        auto proj = gk.dot(grad) * grad;
        //        auto orth = gk - proj;
        //        f_visco += mu * orth;


        //        //if(std::abs(diff.x()) < std::abs(diff.y()))
        //        //    f_visco.x() = f_visco.x() + mu * gk.x();
        //        //else
        //        //    f_visco.y() = f_visco.y() + mu * gk.y();



        //        //f_visco.x() = 0.0;
        //    }
        //}

        fluidVelocity[i] += f_visco;
        //particles[i].vel += 8. * f_visco * (area) / particles[i].rho * dt;
    }
}