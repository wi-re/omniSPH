#include <simulation/SPH.h>
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
#include <random>
#include <iterator>

std::atomic<int64_t> uidCounter = 0;
std::atomic<int64_t> emitCounter = 0;


std::vector<int32_t>& getCell(scalar x, scalar y) {
    if (x != x)
        x = (scalar)0.0;
    if (y != y)
        y = (scalar)0.0;
    std::size_t xi = static_cast<std::size_t>(::floor(std::clamp(x, (scalar)0.0, domainWidth) / scale));
    std::size_t yi = static_cast<std::size_t>(::floor(std::clamp(y, (scalar)0.0, domainHeight) / scale));
    xi = std::clamp(xi, (std::size_t)0, cellsX - 1);
    yi = std::clamp(yi, (std::size_t)0, cellsY - 1);
    return cellArray[yi * (cellsX)+xi];
}

std::pair<int32_t, int32_t> getCellIdx(scalar x, scalar y) {
    if (x != x)
        x = (scalar)0.0;
    if (y != y)
        y = (scalar)0.0;
    std::size_t xi = static_cast<std::size_t>(::floor(std::clamp(x, (scalar)0.0, domainWidth) / scale));
    std::size_t yi = static_cast<std::size_t>(::floor(std::clamp(y, (scalar)0.0, domainHeight) / scale));
    xi = std::clamp(xi, (std::size_t)0, cellsX - 1);
    yi = std::clamp(yi, (std::size_t)0, cellsY - 1);
    return std::make_pair((int32_t)xi, (int32_t)yi);
}

std::vector<int32_t>& getCell(int32_t xi, int32_t yi) {
    return cellArray[yi * (cellsX)+xi];
}

void fillCells() {
    for (int32_t i = 0; i < particles.size(); ++i)
        getCell(particles[i].pos.x(), particles[i].pos.y()).push_back(i);
}

void neighborList() {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());

        for (int32_t xi = -1; xi <= 1; ++xi) {
            for (int32_t yi = -1; yi <= 1; ++yi) {
                const auto& cell = getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    auto& pj = particles[j];
                    vec r = pj.pos - p.pos;
                    if (r.squaredNorm() <= support * support)
                        p.neighbors.push_back(j);
                }
            }
        }
    }
}

void resetFrame() {

    for (auto& v : cellArray)
        v.clear();
    for (auto& p : particles)
        p.reset();

    if (outletSwitch) {
        std::vector<Particle> filteredParticles;
        filteredParticles.reserve(particles.size());
        static auto& t = ParameterManager::instance().get<scalar>("sim.time");
        for (auto& p : particles) {
            if (p.pos.x() > 0.0 && p.pos.x() < domainWidth - 1.25 * domainEpsilon &&
                p.pos.y() > 0.0 && p.pos.y() < domainHeight)
                filteredParticles.push_back(p);
            else {
                if (t == 0.0) continue;
                auto ec = --emitCounter;
                if (ec <= 0) {
                    emitCounter++;
                    filteredParticles.push_back(p);
                }
            }
        }
        particles = filteredParticles;
    }
    std::vector<Particle> filteredParticles;
    filteredParticles.reserve(particles.size());
    static auto& t = ParameterManager::instance().get<scalar>("sim.time");
    for (auto& p : particles) {
        if (p.pos.x() < 0.9 * domainEpsilon || p.pos.x() > domainWidth - 0.9 * domainEpsilon ||
            p.pos.y() < 0.9 * domainEpsilon || p.pos.y() > domainHeight - 0.9 * domainEpsilon);
        else {
            filteredParticles.push_back(p);
        }
    }
    particles = filteredParticles;

    if (particlesDFSPH.size() != particles.size())
        particlesDFSPH.resize(particles.size());
    for (auto& d : particlesDFSPH)
        d.reset();
    static auto& numPtcls = ParameterManager::instance().get<std::size_t>("props.numptcls");
    numPtcls = particles.size();
}

void emitParticles() {
    if (!inletSwitch)
        return;
    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    //auto dt = ParameterManager::instance().get<scalar>("sim.dt");
    static scalar offset = -dt * speed;
    offset = std::fmod(offset + dt * speed, 2.0 * packing_2D);

    auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
    auto t = 12.5 / domainScale;
    auto eps = scale / sqrt(2) * 0.2;
    auto particlesGen = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainEpsilon + 5.0 / domainScale, domainHeight - domainEpsilon - spacing_2D));
    //auto particlesGen = genParticles(vec(15.5 / domainScale, m - t), vec(19.5 / domainScale, m + t));
    for (auto& p : particlesGen) {
        auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());
        p.pos.x() += offset;
        bool emit = true;
        //std::cout << p.pos.x() << " x " << p.pos.y() << " -> " << ix << " : " << iy << std::endl;
        for (int32_t xi = -1; xi <= 1; ++xi) {
            for (int32_t yi = -1; yi <= 1; ++yi) {
                const auto& cell = getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    auto& pj = particles[j];
                    vec r = pj.pos - p.pos;
                    if (r.squaredNorm() <= 0.5 * 2.0 * 2.0 * packing_2D * packing_2D) {
                        emit = false;
                        //pj.vel = vec(speed, 0.0);
                        //pj.angularVelocity = 0.0;
                    }
                }
            }
        }


        p.vel = vec(speed, 0);
        if (emit) {
            p.uid = uidCounter++;
            emitCounter++;
            particles.push_back(p);
            particlesDFSPH.push_back(dfsphState{});
        }
    }
    //std::cout << emitCounter << " [ " << uidCounter << " ] " << std::endl;

}
void initializeParameters(int32_t scene) {
    static auto vectorModes = std::vector<detail::iAny>{
            std::string("magnitude"),
            std::string("x"),
            std::string("y") };
    static auto buffers = std::vector<detail::iAny>{
                std::string("velocity"),
                std::string("angularVelocity"),
                std::string("acceleration"),
                std::string("density"),
                std::string("neighbors"),
                std::string("UID"),
                std::string("alpha"),
                std::string("area"),
                std::string("pressure1"),
                std::string("pressure2"),
                std::string("pressureBoundary"),
                std::string("source"),
                std::string("dpdt"),
                std::string("rhoStar"),
                std::string("predictedVelocity"),
                std::string("pressureAcceleration") };

    ParameterManager::instance().newParameter("colorMap.vectorMode", std::string("magnitude"), { .constant = false 
        });
    ParameterManager::instance().newParameter("colorMap.buffer", std::string("pressure1"), { .constant = false 
        });



    ParameterManager::instance().newParameter("ray.origin", vec(50, 25), { .constant = false });
    ParameterManager::instance().newParameter("ray.target", vec(5, 5), { .constant = false });
    ParameterManager::instance().newParameter("ray.fov", 10.0, { .constant = false , .range = Range{0.1,360.0} });
    ParameterManager::instance().newParameter("ray.resolution", 8, { .constant = false , .range = Range{16,4096} });
    ParameterManager::instance().newParameter("ray.render", false, { .constant = false });
    ParameterManager::instance().newParameter("ray.renderImplicit", true, { .constant = false });
    ParameterManager::instance().newParameter("ray.renderExplicit", false, { .constant = false });
    ParameterManager::instance().newParameter("ray.hScale", 1.0, { .constant = false });
    ParameterManager::instance().newParameter("ray.subSteps", 32, { .constant = false, .range = Range{2,64} });
    ParameterManager::instance().newParameter("ray.gScaleFine", 32.0, { .constant = false, .range = Range{0.01,2.0} });
    ParameterManager::instance().newParameter("ray.gScaleExplicit", 1.0, { .constant = false, .range = Range{0.01,2.0} });
    ParameterManager::instance().newParameter("ray.gScaleImplicit", 1.0, { .constant = false, .range = Range{0.01,2.0} });
    ParameterManager::instance().newParameter("ray.iso", -scalar(0.5), { .constant = false, .range = Range{-1.0, 1.0} });

    ParameterManager::instance().newParameter("sim.frame", 0, { .constant = false });
    ParameterManager::instance().newParameter("field.render", false, { .constant = false });
    ParameterManager::instance().newParameter("field.min", 0.0, { .constant = false , .range = Range{-10.0,10.0} });
    ParameterManager::instance().newParameter("field.max", 1.0, { .constant = false, .range = Range{-10.0,10.0} });
    ParameterManager::instance().newParameter("marching.render", false, { .constant = false });
    ParameterManager::instance().newParameter("marching.solid", false, { .constant = false });
    ParameterManager::instance().newParameter("marching.iso", -scalar(0.5), { .constant = false, .range = Range{-1.0, 1.0} });
    ParameterManager::instance().newParameter("marching.gScale", scalar(1), { .constant = false, .range = Range{0.01,2.0} });
    ParameterManager::instance().newParameter("marching.hScale", scalar(1), { .constant = false, .range = Range{0.01,2.0} });
    ParameterManager::instance().newParameter("marching.method", int32_t(0), { .constant = false, .range = Range{0,0} });
    ParameterManager::instance().newParameter("render.showGrid", false, { .constant = false });

    ParameterManager::instance().newParameter("colorMap.min", scalar(-1.5), { .constant = false , .range = Range{-10.0,10.0} });
    ParameterManager::instance().newParameter("colorMap.map", 0, { .constant = false , .range = Range{0,3} });
    ParameterManager::instance().newParameter("colorMap.limit", false, { .constant = false });
    ParameterManager::instance().newParameter("colorMap.auto", true, { .constant = false });

    ParameterManager::instance().newParameter("colorMap.max", scalar(1.5), { .constant = false , .range = Range{-10.0,10.0} });

    ParameterManager::instance().newParameter("vorticity.nu_t", 0.025, { .constant = false, .range = Range{0.0, 1.0} });
    ParameterManager::instance().newParameter("vorticity.inverseInertia", 0.5, { .constant = false, .range = Range{0.0, 1.0} });
    ParameterManager::instance().newParameter("vorticity.angularViscosity", 0.005, { .constant = false, .range = Range{0.0, 1.0} });

    ParameterManager::instance().newParameter("sim.inletSpeed", 0.0, { .constant = false, .range = Range{0.01,2.00} });
    ParameterManager::instance().newParameter("sim.inletSpeedGoal", .75, { .constant = false, .range = Range{0.01,2.00} });
    ParameterManager::instance().newParameter("sim.time", simulationTime, { .constant = true });
    ParameterManager::instance().newParameter("props.numptcls", particles.size(), { .constant = true });
    ParameterManager::instance().newParameter("props.scale", scale, { .constant = true });
    ParameterManager::instance().newParameter("ptcl.render", true, { .constant = false });
    ParameterManager::instance().newParameter("ptcl.radius", radius, { .constant = true });
    ParameterManager::instance().newParameter("sim.minDt", 0.00001, { .constant = false, .range = Range{0.0001, 0.016} });
    ParameterManager::instance().newParameter("sim.maxDt", 0.002, { .constant = false, .range = Range{0.0001, 0.016} });
    ParameterManager::instance().newParameter("sim.dt", dt, { .constant = true });
    ParameterManager::instance().newParameter("props.domainWidth", domainWidth, { .constant = true });
    ParameterManager::instance().newParameter("props.domainHeight", domainHeight, { .constant = true });
    ParameterManager::instance().newParameter("props.cellsX", cellsX, { .constant = true });
    ParameterManager::instance().newParameter("props.cellsY", cellsY, { .constant = true });
    ParameterManager::instance().newParameter("props.packing_2D", packing_2D, { .constant = true });
    ParameterManager::instance().newParameter("props.support", support, { .constant = true });
    ParameterManager::instance().newParameter("props.damping", 0.0, { .constant = false, .range = Range{0.00,1.00} });
    ParameterManager::instance().newParameter("ptcl.area", area, { .constant = true });
    ParameterManager::instance().newParameter("ptcl.mass", mass, { .constant = true });
    ParameterManager::instance().newParameter("ptcl.maxVelocity", scalar(0), { .constant = true });
    ParameterManager::instance().newParameter("ptcl.viscosityConstant", 0.0001, { .constant = false, .range = Range{0.01, 0.05} });

    ParameterManager::instance().newParameter("dfsph.densityEta", scalar(0.001), { .constant = true });
    ParameterManager::instance().newParameter("dfsph.divergenceEta", scalar(0.001), { .constant = true });
    ParameterManager::instance().newParameter("dfsph.densityIterations", 0, { .constant = true });
    ParameterManager::instance().newParameter("dfsph.divergenceIterations", 0, { .constant = true });
    ParameterManager::instance().newParameter("dfsph.densityError", scalar(0), { .constant = true });
    ParameterManager::instance().newParameter("dfsph.divergenceError", scalar(0), { .constant = true });
}

void initializeSPH(int32_t scene) {
    triangles.clear();
    //gravity = vec(0.0,0.0);
    float d = 1.5;

    //switch (simulationCase) {
    //case cornerAngle::acute:
    //    triangles.push_back(Triangle{ {50, 45}, {50, 25}, {25, 45} }); // 135
    //case cornerAngle::ortho:
    //    triangles.push_back(Triangle{ {50, 25}, {50, 45}, {75, 45} });   // 90
    //case cornerAngle::obtuse:
    //    triangles.push_back(Triangle{ {50, 25}, {95, 25}, {95, 61} });  // 45
    //    triangles.push_back(Triangle{ {50, 25}, {75, 25}, {75, 5} }); // 0
    //case cornerAngle::nobtuse:
    //    triangles.push_back(Triangle{ {50, 25}, {75, 5}, {50, 5} }); // -45
    //case cornerAngle::northo:
    //    triangles.push_back(Triangle{ {50, 25}, {50, 5}, {25, 5} }); // -90
    //case cornerAngle::nacute:
    //    triangles.push_back(Triangle{ {50, 25}, {5, -11}, {5, 25} }); // -135
    //    break;
    //case cornerAngle::Box1:
    //    d = 0.02559;
    //    break;
    //case cornerAngle::Box4:
    //    d = 4.;
    //    break;
    //case cornerAngle::Box1_4:
    //    d = 0.25;
    //    break;
    //default: break;
    //}
    //if (simulationCase == cornerAngle::Box1 || cornerAngle::Box4 == simulationCase || cornerAngle::Box1_4 == simulationCase) {
    //    triangles.push_back(Triangle{ {50,5},{50,5.0 + d},{50.0 + d,5.0 + d} });
    //    triangles.push_back(Triangle{ {50,5},{50.0 + d,5.0 + d},{50.0 + d,5} });

    //}


    //if (simulationCase != cornerAngle::Box1 && simulationCase != cornerAngle::Box1_4 && simulationCase != cornerAngle::Box4) {
    //    vec P1(5, 25);
    //    vec P2(50, 25);
    //    vec P3(0, 0);
    //    if (simulationCase == cornerAngle::acute)
    //        P3 = vec(0, 75);
    //    if (simulationCase == cornerAngle::ortho)
    //        P3 = vec(50, 75);
    //    if (simulationCase == cornerAngle::obtuse)
    //        P3 = vec(100, 75);
    //    if (simulationCase == cornerAngle::nobtuse)
    //        P3 = vec(100, -75);
    //    if (simulationCase == cornerAngle::northo)
    //        P3 = vec(50, -75);
    //    if (simulationCase == cornerAngle::nacute)
    //        P3 = vec(0, -75);

    //    std::vector<vec> polygon;
    //    polygon.push_back(P1);
    //    polygon.push_back(P2);
    //    if (simulationCase == cornerAngle::acute) {
    //        polygon.push_back(vec(25, 45));
    //        polygon.push_back(vec(5, 45));
    //    }
    //    if (simulationCase == cornerAngle::ortho) {
    //        polygon.push_back(vec(50, 45));
    //        polygon.push_back(vec(5, 45));
    //    }
    //    if (simulationCase == cornerAngle::obtuse) {
    //        polygon.push_back(vec(75, 45));
    //        polygon.push_back(vec(5, 45));
    //    }
    //    if (simulationCase == cornerAngle::nobtuse) {
    //        polygon.push_back(vec(75, 5));
    //        polygon.push_back(vec(95, 5));
    //        polygon.push_back(vec(95, 45));
    //        polygon.push_back(vec(5, 45));
    //    }
    //    if (simulationCase == cornerAngle::northo) {
    //        polygon.push_back(vec(50, 5));
    //        polygon.push_back(vec(95, 5));
    //        polygon.push_back(vec(95, 45));
    //        polygon.push_back(vec(5, 45));
    //    }
    //    if (simulationCase == cornerAngle::nacute) {
    //        polygon.push_back(vec(25, 5));
    //        polygon.push_back(vec(95, 5));
    //        polygon.push_back(vec(95, 45));
    //        polygon.push_back(vec(5, 45));
    //    }
    //}
    //else {
    //    vec P1(5, 5);
    //    vec P2(95, 5);
    //    vec P3(95, 45);
    //    vec P4(5, 45);

    //    polygon.push_back(P2);
    //    polygon.push_back(P3);
    //    polygon.push_back(P4);
    //    polygon.push_back(P1);
    //    float d = 0.0559f;
    //    if (simulationCase == cornerAngle::Box4)
    //        d = 4.f;
    //    if (simulationCase == cornerAngle::Box1_4)
    //        d = 0.25f;
    //    if (simulationMethod == boundaryMethod::sdf) {
    //        polygon.push_back(vec{ 50,5 });
    //        polygon.push_back(vec{ 50,5 + d });
    //        polygon.push_back(vec{ 50 + d,5 + d });
    //        polygon.push_back(vec{ 50 + d,5 });
    //    }
    //}
    // thin domain
    vec P1(domainEpsilon, domainEpsilon);
    vec P2(domainWidth - domainEpsilon, domainEpsilon);
    vec P3(domainWidth - domainEpsilon, domainHeight - domainEpsilon);
    vec P4(domainEpsilon, domainHeight - domainEpsilon);

    polygon.push_back(P2);
    polygon.push_back(P3);
    polygon.push_back(P4);
    polygon.push_back(P1);

    //vec P1(domainEpsilon, domainEpsilon);
    //vec P2(domainWidth / 2.0, domainEpsilon);
    //vec P3(domainWidth - domainEpsilon, domainHeight/5.0);
    //vec P4(domainWidth - domainEpsilon, domainHeight - domainEpsilon);
    //vec P5(domainEpsilon, domainHeight - domainEpsilon);
    //polygon.push_back(P2);
    //polygon.push_back(P3);
    //polygon.push_back(P4);
    //polygon.push_back(P5);
    //polygon.push_back(P1);
    //float d = 0.0559f;

    //if (simulationCase == cornerAngle::Box4)
    //    d = 4.f;
    //if (simulationCase == cornerAngle::Box1_4)
    //    d = 0.25f;
    //if (simulationMethod == boundaryMethod::sdf) {
    //    polygon.push_back(vec{ 50,5 });
    //    polygon.push_back(vec{ 50,5 + d });
    //    polygon.push_back(vec{ 50 + d,5 + d });
    //    polygon.push_back(vec{ 50 + d,5 });
    //}

    auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
    auto t = 12.5 / domainScale;

    vec begin(domainEpsilon + 10 / domainScale, m - t);
    vec end(domainWidth - 20 / domainScale, m + t);
    vec mid = (end + begin) * 0.5;
    mid.x() -= 25.0 / domainScale;
    mid.y() -= 0.1 / domainScale;

    auto thickness = 1.0 / domainScale;
    auto c = 3.5 / domainScale;


    std::vector<vec> Upper{ {begin.x(), begin.y()},
                           {end.x(), begin.y()},
                           {end.x(), begin.y() - thickness},
                           {begin.x(), begin.y() - thickness} };
    std::vector<vec> Lower{ {begin.x(), end.y()},
                           {end.x(), end.y()},
                           {end.x(), end.y() + thickness},
                           {begin.x(), end.y() + thickness} };
    std::vector<vec> Obs{ {mid.x() - c, mid.y() - 0.75 * c},
                           {mid.x() + 4.0 * c, mid.y() - 2.0 * c},
                           {mid.x() + 4.0 * c, mid.y() + 2.0 * c},
                           {mid.x() - c, mid.y() + 0.75 * c} };
    std::vector<vec> ObsBoxCenter{ {mid.x() - c, mid.y() - c},
                           {mid.x() + c, mid.y() - c},
                           {mid.x() + c, mid.y() + c},
                           {mid.x() - c, mid.y() + c} };

    std::vector<vec> ObsBox{ {mid.x() - c, domainEpsilon},
                           {mid.x() + c, domainEpsilon},
                           {mid.x() + c, domainEpsilon + c},
                           {mid.x() - c, domainEpsilon + c} };

    //triangles.push_back(Triangle{ {begin.x(),begin.y()},{end.x(),begin.y()},{end.x(),begin.y()- thickness} });
    //triangles.push_back(Triangle{ {begin.x(),begin.y()},{end.x(),begin.y()- thickness},{begin.x(),begin.y() - thickness} });
    //triangles.push_back(Triangle{ {begin.x(),end.y()},{end.x(),end.y()},{end.x(),end.y() + thickness} });
    //triangles.push_back(Triangle{ {begin.x(),end.y()},{end.x(),end.y() + thickness},{begin.x(),end.y() + thickness} }); 
    //triangles.push_back(Triangle{{mid.x() - c, mid.y() - c},
    //                             {mid.x() - c, mid.y() + c},
    //                             {mid.x() + c, mid.y() + c}});
    //triangles.push_back(Triangle{{mid.x() - c, mid.y() - c},
    //                             {mid.x() + c, mid.y() + c},
    //                             {mid.x() + c, mid.y() - c}}); 
    auto moffset = 60.0 / domainScale;
    vec mm = (end + begin) * 0.5;
    auto width = 25.0 / domainScale;
    auto height = 10.0 / domainScale;

    auto rows = height / packing_2D;
    height = ::ceil(rows) * packing_2D;

    mm.x() -= moffset;


    auto left = mm.x() - domainEpsilon - width - 2. * spacing_2D;
    auto columns = left / packing_2D;

    std::cout << "Packing: " << packing_2D << std::endl;
    std::cout << "Obstacle Center: " << mm.x() << std::endl;
    std::cout << "Obstacle Left Edge: " << mm.x() - width << std::endl;
    std::cout << "Left Gap: " << mm.x() - domainEpsilon - width - 2. * spacing_2D << " [ " << left << " ]" << std::endl;
    std::cout << "Left Columns: " << columns << std::endl;
    std::cout << "New Columns: " << ::ceil(columns) << std::endl;
    left = domainEpsilon + ::ceil(columns) * packing_2D + 2. * spacing_2D;
    width = ::ceil(2.0 * width / packing_2D) * packing_2D;
    std::cout << "New Gap: " << left << " [ " << ::ceil(left) * packing_2D << " + " << domainEpsilon << " + " << 2 * spacing_2D << " ]" << std::endl;
    std::cout << "New Width: " << width << " was " << 25.0 / domainScale << std::endl;

    auto centerOffset = mm.x();

    std::vector<vec> ObsB{ {left, domainEpsilon},
                           {left + width, domainEpsilon},
                           {left + width, domainEpsilon + height},
                           {left, domainEpsilon + height} };

    std::vector<vec> ObsT{ {left, domainHeight - domainEpsilon},
                           {left + width, domainHeight - domainEpsilon},
                           {left + width, domainHeight - domainEpsilon - height},
                           {left, domainHeight - domainEpsilon - height} };


    std::vector<vec> ObsB2{ {left, domainEpsilon},
                           {left + width, domainEpsilon},
                           {left + width, domainEpsilon + 1.5 * height - 5.0 / domainScale},
                           {left, domainEpsilon + height - 5.0 / domainScale} };

    std::vector<vec> ObsT2{ {left, domainHeight - domainEpsilon},
                           {left + width, domainHeight - domainEpsilon},
                           {left + width, domainHeight - domainEpsilon - 1.5 * height},
                           {left, domainHeight - domainEpsilon - height} };

    std::vector<vec> ObsStair{ {domainEpsilon,domainHeight/5.0},
                           {domainEpsilon + 75.0 / domainScale, domainHeight / 5.0},
                           {domainEpsilon + 75.0 / domainScale, domainHeight / 5.0 - 2.5 * scale},
                           {domainEpsilon, domainHeight / 5.0 - 2.5 * scale} };
    //obstacles.push_back(ObsStair);

    //obstacles.push_back(Upper);
    //obstacles.push_back(Lower);
    //obstacles.push_back(Obs);
    // 
    // thin domain obstacles
    //switch (bc) {
    //case boundaryConfig::Box:
    //    obstacles.push_back(ObsBox); break;
    //case boundaryConfig::CenterBox:
    //    obstacles.push_back(ObsBoxCenter); break;
    //case boundaryConfig::ObstacleBT:
    //    obstacles.push_back(ObsB);
    //    obstacles.push_back(ObsT); break;
    //case boundaryConfig::ObstacleBTLow:
    //    obstacles.push_back(ObsB2);
    //    obstacles.push_back(ObsT2); break;
    //case boundaryConfig::Trapezoid:
    //    obstacles.push_back(Obs); break;
    //}
    //obstacles.push_back(ObsM);
    //obstacles.push_back(ObsB);
    //obstacles.push_back(ObsT);

//	std::cout << boundaryKernel(-1.0) << " : " << boundaryKernel(0.0) << " : " << boundaryKernel(1.0) << std::endl;
  // particles = genParticles(vec(domainEpsilon+0.5, domainEpsilon+0.5), vec(domainWidth / 8, domainHeight -
  // domainEpsilon * 2.0));
  //auto particles1 = genParticles(vec(27.5, 7.5), vec(47.5, 27.5));
   // auto particles1 = genParticles(vec(15, 35), vec(15.1, 35.1));
   ////auto particles1 = genParticles(vec(20, 20), vec(30, 30));
   // //auto particles2 = genParticles(vec(10.25, 10.25), vec(20.25, 20.25));
   // //auto particles2 = genParticles(vec(10.25, 10.25), vec(10.5, 10.5));
   // auto particles2 = genParticles(vec(5 + offset, 5 + offset), vec(20 + offset, 25 + offset));
   // auto particles4 = genParticles(vec(45.25, 25.25), vec(55.25, 35.25));
   // auto particles3 = genParticles(vec(5.25 / domainScale, 25.25), vec(20, 40));
   // auto offset = domainEpsilon / 5;

   // auto eps = scale / sqrt(2) * 0.2;


    if (pc != particleConfig::None) {
        std::vector<Particle> particles5;
        // thin domain
        //switch (pc) {
        //case particleConfig::Domain:particles5 = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainWidth - domainEpsilon - spacing_2D, domainHeight - domainEpsilon - spacing_2D)); break;
        //case particleConfig::DamBreak:particles5 = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainEpsilon + 20.0 / domainScale, domainEpsilon + 37.5 / domainScale)); break;
        //}
        particles5 = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainEpsilon + 40.0 / domainScale, domainHeight * 3.0 / 5.0));
        //auto candidates = particles3;
        //if (simulationCase != cornerAngle::Box1 && simulationCase != cornerAngle::Box1_4 && simulationCase != cornerAngle::Box4) {
        //    candidates = particles3;
        //}
        //else {
        //    candidates = particles2;

        //}
        auto candidates = particles5;

        for (int32_t i = 0; i < candidates.size(); ++i)
            getCell(candidates[i].pos.x(), candidates[i].pos.y()).push_back(i);

        auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
        std::vector<Particle> toEmit;
        for (auto& p : candidates) {
            auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());
            bool emit = true;
            auto density = 0.0;
            for (int32_t xi = -1; xi <= 1; ++xi) {
                for (int32_t yi = -1; yi <= 1; ++yi) {
                    const auto& cell = getCell(ix + xi, iy + yi);
                    for (auto j : cell) {
                        auto& pj = candidates[j];
                        vec r = pj.pos - p.pos;
                        density += area * W(p.pos, pj.pos);
                        //pi.rho = std::max(pi.rho, 0.5);
                    }
                }
            }
            boundaryFunc(p.pos, [&density](auto bpos, auto d, auto k, auto gk, auto triangle) {
                //std::cout << i << " - " << pi.pos.x() << " : " << pi.pos.y() << " -> " << bpos.x() << " : " << bpos.y() << ", " << d << ", " << k << ", " << gk.x() << " : " << gk.y() << std::endl;
                if (d >= -1.1 * spacing_2D)
                    density += k;
                });
            if (density < 1.02) {
                p.uid = uidCounter++;
                auto dist = p.pos.x() - 2.0 * domainEpsilon;
                //if (dist < 35.0 / domainScale) {
                //    dist /= 35.0 / domainScale;
                //    p.vel = vec((1.0-dist) * speed, 0.0);
                //}

                toEmit.push_back(p);
            }
        }

        std::default_random_engine engineObj(std::chrono::system_clock::now().time_since_epoch().count());
        for (int32_t i = 0; i < 250; ++i) {
            //creat a random engine object and seed it 
            //The distribution will be used
            std::uniform_int_distribution<size_t> randomInt{ 0ul, toEmit.size() - 1 };

            //Erase
            //toEmit.erase(toEmit.begin() + randomInt(engineObj));
        }
        for (auto& p : toEmit)
            particles.push_back(p);

    }
    //if (simulationCase != cornerAngle::Box1 && simulationCase != cornerAngle::Box1_4 && simulationCase != cornerAngle::Box4) {
    //    for (auto& p : particles3)
    //        particles.push_back(p);
    //}
    //else {
    //    for (auto& p : particles2) {
    //        //p.vel = vec(5, 0);
    //        particles.push_back(p);
    //    }
    //    for (auto& p : particles4) {
    //        p.vel = vec(-5, 0);
    //       // particles.push_back(p);
    //    }

    //}
    //triangles.push_back(Triangle{ {5,5},{95,5},{95,0} });
    //triangles.push_back(Triangle{ {5,5},{95,0},{5,0} });
    //triangles.push_back(Triangle{ {0,0},{0,50},{5,50} });
    //triangles.push_back(Triangle{ {0,0},{5,50},{5,0} });



}