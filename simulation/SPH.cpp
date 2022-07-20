
//#define _CRT_SECURE_NO_WARNINGS
#include "SPH.h"
#include "2DMath.h"
#include <algorithm>
#include <array>
#include <boost/range/combine.hpp>
#include <iostream>
#include <numeric>
#include <chrono>
#include <sstream>
#include <atomic>


#include <time.h>
#include <iomanip>
#include <tools/timer.h>
#include <simulation/2DMath.h>
#include <cfloat>

using clk = std::chrono::high_resolution_clock;
scalar toMs(clk::duration dur) {
  return static_cast<scalar>(std::chrono::duration_cast<std::chrono::microseconds>(dur).count()) / scalar(1000.0);
}



void SPHSimulation::timestep(){
    TIME_CODE(0, "Simulation - Overall",
        TIME_CODE(1, "Simulation - Reset", resetFrame());
        TIME_CODE(2, "Simulation - Cell construction", fillCells());
        TIME_CODE(3, "Simulation - Neighbor search", neighborList());
        TIME_CODE(4, "Simulation - Density", density());
        TIME_CODE(5, "Simulation - Vorticity", computeVorticity());
        TIME_CODE(6, "Simulation - External", externalForces());
        // if(!useEOS){
             TIME_CODE(7, "Simulation - Divergence", divergenceSolve());
             TIME_CODE(8, "Simulation - Density", densitySolve());
        // }
        // else{
        //     TIME_CODE(7, "Simulation - EOS Pressure", stateSPH());
        // }
         TIME_CODE(9, "Simulation - XSPH", XSPH());
        //TIME_CODE(10, "Simulation - Vorticity", refineVorticity());
        TIME_CODE(11, "Simulation - Integration", Integrate());
        //TIME_CODE(12, "Simulation - BXSPH", BXSPH());
        TIME_CODE(13, "Simulation - Dump", dump());
        TIME_CODE(14, "Simulation - Emission", emitParticles());
    );
}

void SPHSimulation::fillCells(){    
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");

    for (int32_t i = 0; i < numPtcls; ++i)
        getCell(fluidPosition[i].x(), fluidPosition[i].y()).push_back(i);
}
void SPHSimulation::neighborList(){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        auto& p = fluidPosition[i];
        auto& neighs = fluidNeighborList[i];
        auto& hi = fluidSupport[i];
        auto [ix, iy] = getCellIdx(p.x(), p.y());

        for (int32_t xi = -1; xi <= 1; ++xi) {
            for (int32_t yi = -1; yi <= 1; ++yi) {
                auto xxi = std::clamp(ix + xi, 0, (int32_t) cellsX -1);
                auto yyi = std::clamp(iy + yi, 0, (int32_t) cellsY -1);
                const auto& cell = getCell(xxi, yyi);
                for (auto j : cell) {
                    auto& pj = fluidPosition[j];
                    vec r = pj - p;
                    auto hj = fluidSupport[j];
                    auto support = (hi + hj) / 2.;
                    if (r.squaredNorm() <= support * support)
                        neighs.push_back(j);
                }
            }
        }
    }

    
    auto triangleArea = [](Triangle t) {
        auto [p0, p1, p2] = t;
        return 0.5 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1]) + epsilon;
    };
    auto pointInTriangle = [triangleArea](vec p, Triangle tri) {
        auto [p0, p1, p2] = tri;
        auto area = triangleArea(tri);
        auto s = 1.0 / (2.0 * area) * (p0[1] * p2[0] - p0[0] * p2[1] + (p2[1] - p0[1]) * p[0] + (p0[0] - p2[0]) * p[1]);
        auto t = 1.0 / (2.0 * area) * (p0[0] * p1[1] - p0[1] * p1[0] + (p0[1] - p1[1]) * p[0] + (p1[0] - p0[0]) * p[1]);
        auto u = 1.0 - s - t;
        if (s > 0.0 && t > 0.0 && u > 0.0) {
            return 1.0;
        }
        if (s >= 0.0 && t >= 0.0 && u >= 0.0) {
            return 0.0;
        }
        return -1.0;
    };
    auto pointInCircle = [](vec c, scalar r, vec p) {
        auto a = p.x() - c.x();
        auto b = p.y() - c.y();
        auto ab2 = a * a + b * b;
        auto r2 = r * r;
        if (ab2 < r2)
            return 1.0;
        if (abs(ab2) / r2 < 1e-11)
            return 0.0;
        return -1.0;
    };
    auto closestPoint = [](vec P, vec A, vec B, bool clipped) {
        vec ap = P - A;
        vec ab = B - A;
        scalar ab2 = ab.dot(ab) + epsilon;
        scalar apab = ap.dot(ab);
        scalar t = apab / ab2;
        if (clipped)
            t = std::clamp(t, 0.0, 1.0);
        return (vec)(A + ab * t);
    };
    auto closestPointEdge = [closestPoint](vec c, vec p1, vec p2, vec center, bool check) {
        auto dC = (p2[1] - p1[1]) * c[0] - (p2[0] - p1[0]) * c[1] + p2[0] * p1[1] - p2[1] * p1[0];
        auto dT = (p2[1] - p1[1]) * center[0] - (p2[0] - p1[0]) * center[1] + p2[0] * p1[1] - p2[1] * p1[0];
        if (dC * dT < 0.0 || !check)
            return closestPoint(c, p1, p2, false);
        return c;
    };
    auto closestPointTriangle = [pointInTriangle, closestPointEdge](vec P, Triangle tri) {
        if (pointInTriangle(P, tri) >= 0) {
            auto triCenter0 = (tri.v0 + tri.v1 + tri.v2) / 3.0;

            auto P01 = closestPointEdge(P, tri.v0, tri.v1, triCenter0, false);
            auto P12 = closestPointEdge(P, tri.v1, tri.v2, triCenter0, false);
            auto P20 = closestPointEdge(P, tri.v2, tri.v0, triCenter0, false);

            auto d01 = (P01 - P).norm();
            auto d12 = (P12 - P).norm();
            auto d20 = (P20 - P).norm();
            if (d01 <= d12 && d01 <= d20)
                return std::make_pair(P01, d01);
            if (d12 <= d01 && d12 <= d20)
                return std::make_pair(P12, d12);
            return std::make_pair(P20, d20);
        }
        else {
            auto triCenter0 = (tri.v0 + tri.v1 + tri.v2) / 3.0;

            auto P01 = closestPointEdge(P, tri.v0, tri.v1, triCenter0, true);
            auto P12 = closestPointEdge(P01, tri.v1, tri.v2, triCenter0, true);
            auto P20 = closestPointEdge(P12, tri.v2, tri.v0, triCenter0, true);

            auto d01 = (P01 - P).norm();
            auto d12 = (P12 - P).norm();
            auto d20 = (P20 - P).norm();
            return std::make_pair(P20, -d20);
            if (d01 <= d12 && d01 <= d20)
                return std::make_pair(P01, d01);
            if (d12 <= d01 && d12 <= d20)
                return std::make_pair(P12, d12);
            return std::make_pair(P20, d20);
        }
    };

#pragma omp parallel for
    for (int32_t i = 0; i < numPtcls; ++i) {
        auto& p = fluidPosition[i];
        auto& support = fluidSupport[i];
        int32_t ti = 0;
        auto [xi, yi] = getCellIdx(p.x(), p.y());
        auto triangleList = getTriangleCell(xi, yi);
            for (const auto& ti : triangleList) {
                const auto& t = boundaryTriangles[ti];
                auto [cP, d] = closestPointTriangle(p, t);
                if (d > -support) {
                    fluidTriangleNeighborList[i].push_back(ti);
                }
            }
    }



}

std::tuple<scalar, scalar, std::vector<scalar>> SPHSimulation::colorMap(property_t prop, bool autoMinMax, scalar min, scalar max){
    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    std::string vectorMode = pm.get<std::string>("colorMap.vectorMode");
    int32_t mode = 0;
    if(vectorMode == "magnitude")
        mode = 0;
    if(vectorMode == "x")
        mode = 1;
    if(vectorMode == "y")
        mode = 2;
    auto mapVec = [mode](vec vector){
        switch(mode){
            case 0: return vector.norm();
            case 1: return vector.x();
            case 2: return vector.y();
        }
    };
    std::vector<scalar> mappedData(numPtcls, 0.);
    if(autoMinMax){
        min = DBL_MAX;
        max = - DBL_MAX;
        for(int32_t i = 0; i < numPtcls; ++ i){
            auto s = 0.;
            switch(prop){
                case property_t::position:          s = mapVec(fluidPosition[i]); break;
                case property_t::velocity:          s = mapVec(fluidVelocity[i]); break;
                case property_t::accel:             s = mapVec(fluidAccel[i]); break;
                case property_t::predVelocity:      s = mapVec(fluidPredVelocity[i]); break;
                case property_t::predAccel:         s = mapVec(fluidPredAccel[i]); break;
                case property_t::predPosition:      s = mapVec(fluidPredPosition[i]); break;
                case property_t::density:           s = fluidDensity[i]; break;
                case property_t::vorticity:         s = fluidVorticity[i]; break;
                case property_t::angularVelocity:   s = fluidAngularVelocity[i]; break;
                case property_t::area:              s = fluidArea[i]; break;
                case property_t::restDensity:       s = fluidRestDensity[i]; break;
                case property_t::support:           s = fluidSupport[i]; break;
                case property_t::alpha:             s = fluidAlpha[i]; break;
                case property_t::actualArea:        s = fluidActualArea[i]; break;
                case property_t::pressure1:         s = fluidPressure1[i]; break;
                case property_t::pressure2:         s = fluidPressure2[i]; break;
                case property_t::boundaryPressure:  s = fluidBoundaryPressure[i]; break;
                case property_t::sourceTerm:        s = fluidSourceTerm[i]; break;
                case property_t::dpdt:              s = fluidDpDt[i]; break;
                case property_t::rhoStar:           s = fluidDensityStar[i]; break;
                case property_t::priorPressure:     s = fluidPriorPressure[i]; break;
                case property_t::UID:               s = fluidUID[i]; break;
                case property_t::neighbors:         s = fluidNeighborList[i].size(); break;
            }
            min = std::min(s, min);
            max = std::max(s, max);
        }
    }

    for(int32_t i = 0; i < numPtcls; ++ i){
            auto s = 0.;
            switch(prop){
                case property_t::position:          s = mapVec(fluidPosition[i]); break;
                case property_t::velocity:          s = mapVec(fluidVelocity[i]); break;
                case property_t::accel:             s = mapVec(fluidAccel[i]); break;
                case property_t::predVelocity:      s = mapVec(fluidPredVelocity[i]); break;
                case property_t::predAccel:         s = mapVec(fluidPredAccel[i]); break;
                case property_t::predPosition:      s = mapVec(fluidPredPosition[i]); break;
                case property_t::density:           s = fluidDensity[i]; break;
                case property_t::vorticity:         s = fluidVorticity[i]; break;
                case property_t::angularVelocity:   s = fluidAngularVelocity[i]; break;
                case property_t::area:              s = fluidArea[i]; break;
                case property_t::restDensity:       s = fluidRestDensity[i]; break;
                case property_t::support:           s = fluidSupport[i]; break;
                case property_t::alpha:             s = fluidAlpha[i]; break;
                case property_t::actualArea:        s = fluidActualArea[i]; break;
                case property_t::pressure1:         s = fluidPressure1[i]; break;
                case property_t::pressure2:         s = fluidPressure2[i]; break;
                case property_t::boundaryPressure:  s = fluidBoundaryPressure[i]; break;
                case property_t::sourceTerm:        s = fluidSourceTerm[i]; break;
                case property_t::dpdt:              s = fluidDpDt[i]; break;
                case property_t::rhoStar:           s = fluidDensityStar[i]; break;
                case property_t::priorPressure:     s = fluidPriorPressure[i]; break;
                case property_t::UID:               s = fluidUID[i]; break;
                case property_t::neighbors:         s = fluidNeighborList[i].size(); break;
            }
            mappedData[i] = (s - min) / (max - min);
            if(max == min)
                mappedData[i] = 0.;
    }
    return std::make_tuple(min, max, mappedData);
}

