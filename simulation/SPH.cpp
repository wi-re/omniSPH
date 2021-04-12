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

constexpr bool inletSwitch = false;
constexpr bool outletSwitch = false;
constexpr bool gravitySwitch = true;
constexpr bool backgroundPressureSwitch = false;
enum struct particleConfig {
    Domain, DamBreak, None
};
enum struct boundaryConfig {
    Box, ObstacleBT, ObstacleBTLow, Trapezoid, CenterBox
};
constexpr boundaryConfig bc = boundaryConfig::Box;
constexpr particleConfig pc = particleConfig::DamBreak;

std::atomic<int64_t> uidCounter = 0;
std::atomic<int64_t> emitCounter = 0;

std::vector<int32_t> &getCell(scalar x, scalar y) {
  if (x != x)
    x = (scalar)0.0;
  if (y != y)
    y = (scalar)0.0;
  std::size_t xi = static_cast<std::size_t>(::floor(std::clamp(x, (scalar)0.0, domainWidth) / scale) );
  std::size_t yi = static_cast<std::size_t>(::floor(std::clamp(y, (scalar)0.0, domainHeight) / scale) );
  xi = std::clamp(xi, (std::size_t)0, cellsX - 1);
  yi = std::clamp(yi, (std::size_t)0, cellsY - 1);
  return cellArray[yi * (cellsX) + xi];
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
    auto &p = particles[i];
    auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());

    for (int32_t xi = -1; xi <= 1; ++xi) {
      for (int32_t yi = -1; yi <= 1; ++yi) {
        const auto &cell = getCell(ix + xi, iy + yi);
        for (auto j : cell) {
          auto &pj = particles[j];
          vec r = pj.pos - p.pos;
          if (r.squaredNorm() <= support * support)
            p.neighbors.push_back(j);
        }
      }
    }
  }
}

void density() {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    for (int32_t j : pi.neighbors) {
        auto& pj = particles[j];
        pi.rho += area * W(pi.pos, pj.pos);
    }
    //pi.rho = std::max(pi.rho, 0.5);
    boundaryFunc(pi.pos, [&pi, i](auto bpos, auto d, auto k, auto gk, auto triangle) {
		//std::cout << i << " - " << pi.pos.x() << " : " << pi.pos.y() << " -> " << bpos.x() << " : " << bpos.y() << ", " << d << ", " << k << ", " << gk.x() << " : " << gk.y() << std::endl;
		pi.rho += k; });
    
  }
}

void computeVorticity() {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        scalar voriticity = 0.0;

        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto grad = gradW(pi.pos, pj.pos);
            auto vel = pi.vel - pj.vel;
            auto term = area / pj.rho * (vel.x() * grad.y() - vel.y() * grad.x());

            voriticity += term;
        }
        pi.vorticity = voriticity;
    }
}

void refineVorticity() {
    static auto& nu_t = ParameterManager::instance().get<scalar>("vorticity.nu_t");
    static auto& intertiaInverse = ParameterManager::instance().get<scalar>("vorticity.inverseInertia");
    static auto& angularViscosity = ParameterManager::instance().get<scalar>("vorticity.angularViscosity");

    std::vector<scalar> dwdt(particles.size(), 0.0);
    std::vector<vec> dvdt(particles.size(), vec(0.0,0.0));
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        scalar voriticity = 0.0;
        vec vterm(0.0,0.0);

        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto grad = gradW(pi.pos, pj.pos);
            auto vel = pi.vel - pj.vel;
            auto vor = pi.angularVelocity - pj.angularVelocity;
            auto term = area / pj.rho * (vel.x() * grad.y() - vel.y() * grad.x());

            voriticity += term;
            vterm.x() += area / pj.rho * (- vor * grad.y());
            vterm.y() += area / pj.rho * (  vor * grad.x());


        }
        dwdt[i] = nu_t * (voriticity - 2.0 * pi.angularVelocity) * intertiaInverse;
        dvdt[i] = nu_t * vterm;
    }

#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        pi.angularVelocity += dwdt[i] * dt;
        pi.vel += dvdt[i] * dt;
    }

#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        scalar angularTerm = 0.0;
        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto kernel = W(pi.pos, pj.pos);
            auto vor = pj.angularVelocity - pi.angularVelocity;
            auto term = area / pj.rho * vor * kernel;
            angularTerm += term;
        }
        dwdt[i] = angularViscosity * angularTerm;
    }
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        pi.angularVelocity += dwdt[i];
    }

//    static auto& nu_t = ParameterManager::instance().get<scalar>("vorticity.nu_t");
//    static auto& intertiaInverse = ParameterManager::instance().get<scalar>("vorticity.inverseInertia");
//    static auto& angularViscosity = ParameterManager::instance().get<scalar>("vorticity.angularViscosity");
//
//    std::vector<scalar> newVorticity(particles.size(), 0.0), estVorticity(particles.size(), 0.0);
//    std::vector<scalar> dissipation(particles.size(), 0.0);
//    std::vector<scalar> streamFn(particles.size(), 0.0);
//    std::vector<scalar> dwdt(particles.size(), 0.0);
//    std::vector<vec> dvdt(particles.size(), vec(0.0,0.0));
//    // Vorticity through linear field
//#pragma omp parallel for
//    for (int32_t i = 0; i < particles.size(); ++i) {
//        auto& pi = particles[i];
//        scalar voriticity = 0.0;
//        vec vterm(0.0, 0.0);
//
//        for (int32_t j : pi.neighbors) {
//            auto& pj = particles[j];
//            auto grad = gradW(pi.pos, pj.pos);
//            auto vel = pi.vel - pj.vel;
//            auto vor = pi.angularVelocity - pj.angularVelocity;
//            auto term = area / pj.rho * (vel.x() * grad.y() - vel.y() * grad.x());
//
//            voriticity += term;
//        }
//        newVorticity[i] = voriticity;
//    }
//    // Vorticity Equation
//#pragma omp parallel for
//    for (int32_t i = 0; i < particles.size(); ++i) {
//        auto& pi = particles[i];
//        scalar voriticity = 0.0;
//        scalar vterm = 0.0;
//
//        for (int32_t j : pi.neighbors) {
//            if (i == j) continue;
//            auto& pj = particles[j];
//            auto grad = gradW(pi.pos, pj.pos);
//            auto diff = pi.pos - pj.pos;
//            auto vel = pi.vel - pj.vel;
//            auto vor = pi.angularVelocity - pj.angularVelocity;
//            auto term = area / pj.rho * vor * grad.norm() / diff.norm();
//
//            voriticity += area / pj.rho * vor * vel.dot(grad);
//            vterm += -area / pj.rho * vor * 2.0 * grad.norm() / (diff.norm() + 0.01 * support * support);
//        }
//        dwdt[i] = voriticity + 0.2 * vterm;
//        estVorticity[i] = pi.vorticity + dwdt[i] * dt;
//        dissipation[i] = estVorticity[i] - newVorticity[i];
//        //std::cout << i << ": " << pi.vorticity << " -> " << estVorticity[i] << " : " << dissipation[i] << " => " << voriticity << " : " << vterm << " --> " << dwdt[i] << std::endl;
//    }
//    // stream function
//#pragma omp parallel for
//    for (int32_t i = 0; i < particles.size(); ++i) {
//        auto& pi = particles[i];
//        scalar voriticity = 0.0;
//        vec vterm(0.0, 0.0);
//
//        for (int32_t j : pi.neighbors) {
//            if (i == j) continue;
//            auto& pj = particles[j];
//            auto diff = pi.pos - pj.pos;
//
//            voriticity += 1.0 / (4.0 * M_PI) * dissipation[j] * area / diff.norm();
//        }
//        streamFn[i] = voriticity;
//    }
//
//#pragma omp parallel for
//    for (int32_t i = 0; i < particles.size(); ++i) {
//        auto& pi = particles[i];
//        vec vterm = vec(0.0,0.0);
//        for (int32_t j : pi.neighbors) {
//            auto& pj = particles[j];
//            auto grad = gradW(pi.pos, pj.pos);
//            auto vor = streamFn[i] - streamFn[j];
//
//            vterm.x() += area / pj.rho * (- vor * grad.y());
//            vterm.y() += area / pj.rho * (  vor * grad.x());
//        }
//        dvdt[i] = vterm;
//        vterm *= 1.0;
//        pi.angularVelocity = pi.vorticity;
//        pi.vel += vterm;
//        //            vterm.x() += area / pj.rho * (- vor * grad.y());
//        //            vterm.y() += area / pj.rho * (  vor * grad.x());
//    }
//


}


void externalForces() {
    if(gravitySwitch)
  for (auto &p : particles)
    p.accel += gravity;
}

void XSPH() {
    static auto& viscosityConstant = ParameterManager::instance().get<scalar>("ptcl.viscosityConstant");
  std::vector<vec> tempV;
//#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    tempV.push_back(pi.vel);
    for (auto &j : pi.neighbors) {
      auto &pj = particles[j];
      tempV[i] += viscosityConstant * (pi.rho + pj.rho) *  area / (pi.rho + pj.rho) * scalar(2.0) * W(pi.pos, pj.pos) * (pj.vel - pi.vel);
    }
  }
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i)
    particles[i].vel = tempV[i];
}

scalar sdpoly(std::vector<vec> v, vec p) {
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

void Integrate(void) {
    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    auto speed2 = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    auto speedGoal = ParameterManager::instance().get<scalar>("sim.inletSpeedGoal");
    static auto& time = ParameterManager::instance().get<scalar>("sim.time");
    auto delay = 0.13;
    if (time < delay)
        speed2 = 0.0;
    else
        speed2 = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    //else if (time < 1.0 + delay) {
    //    speed2 = speedGoal * (time - delay);
    //}

    //speed2 = speed;

    auto particlesGen2 = genParticles(vec(domainWidth - 2.5 * domainEpsilon, domainEpsilon + spacing_2D), vec(domainWidth - domainEpsilon, domainHeight - domainEpsilon - spacing_2D));
    //auto particlesGen2 = genParticles(vec(domainWidth - 24.5 / domainScale, m - t), vec(domainWidth - 19.5 / domainScale, m + t));
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        if(inletSwitch)
        if (p.pos.x() > domainWidth - 2. * domainEpsilon) {
            p.vel = vec(speed2, 0.0);
            p.angularVelocity = 0.0;
            p.accel = vec(0.0,0.0);
        }
        if(outletSwitch)
        if (p.pos.x() < 2. * domainEpsilon) {
            p.vel = vec(speed, 0.0);
            p.angularVelocity = 0.0;
            p.accel = vec(0.0, 0.0);
        }

        //auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());
        //bool emit = true;
        ////std::cout << p.pos.x() << " x " << p.pos.y() << " -> " << ix << " : " << iy << std::endl;
        //for (int32_t xi = -1; xi <= 1; ++xi) {
        //    for (int32_t yi = -1; yi <= 1; ++yi) {
        //        const auto& cell = getCell(ix + xi, iy + yi);
        //        for (auto j : cell) {
        //            auto& pj = particles[j];
        //            vec r = pj.pos - p.pos;
        //            if (r.squaredNorm() <= 0.9 * 2.0 * 2.0 * packing_2D * packing_2D) {
        //                emit = false;
        //            }
        //        }
        //    }
        //}


        ////p.vel = vec(speed, 0);
        ////if (emit) {
        ////    p.uid = uidCounter++;
        ////    particles.push_back(p);
        ////    particlesDFSPH.push_back(dfsphState{});
        ////}
    }
    //auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    ////auto dt = ParameterManager::instance().get<scalar>("sim.dt");
    //static scalar offset = 0.0;
    //offset = std::fmod(offset + dt * speed, packing_2D);

    //auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
    ////auto t = 12.5 / domainScale;
    //auto eps = scale / sqrt(2) * 0.2;
    //auto particlesGen = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainEpsilon + 5.0 / domainScale, domainHeight - domainEpsilon - spacing_2D));
    ////auto particlesGen = genParticles(vec(15.5 / domainScale, m - t), vec(19.5 / domainScale, m + t));
    //for (auto& p : particlesGen) {
    //    auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());
    //    bool emit = true;
    //    //std::cout << p.pos.x() << " x " << p.pos.y() << " -> " << ix << " : " << iy << std::endl;
    //    for (int32_t xi = -1; xi <= 1; ++xi) {
    //        for (int32_t yi = -1; yi <= 1; ++yi) {
    //            const auto& cell = getCell(ix + xi, iy + yi);
    //            for (auto j : cell) {
    //                auto& pj = particles[j];
    //                vec r = pj.pos - p.pos;
    //                if (r.squaredNorm() <= 0.99 * 2.0 * 2.0 * packing_2D * packing_2D) {
    //                    emit = false;
    //                    pj.vel = vec(speed, 0.0);
    //                    pj.angularVelocity = 0.0;
    //                }
    //            }
    //        }
    //    }


    //    p.vel = vec(speed, 0);
    //}

    //auto speed2 = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    //auto speedGoal = ParameterManager::instance().get<scalar>("sim.inletSpeedGoal");
    //static auto& time = ParameterManager::instance().get<scalar>("sim.time");
    //if (time < 0.25)
    //    speed2 = 0.0;
    //else if (time < 2.0) {
    //    speed2 = speedGoal * time / 2.0;
    //}

    //speed2 = speed;

    //auto particlesGen2 = genParticles(vec(domainWidth - 2.5 * domainEpsilon, domainEpsilon + spacing_2D), vec(domainWidth - domainEpsilon, domainHeight - domainEpsilon - spacing_2D));
    ////auto particlesGen2 = genParticles(vec(domainWidth - 24.5 / domainScale, m - t), vec(domainWidth - 19.5 / domainScale, m + t));
    //for (auto& p : particlesGen2) {
    //    auto [ix, iy] = getCellIdx(p.pos.x(), p.pos.y());
    //    bool emit = true;
    //    //std::cout << p.pos.x() << " x " << p.pos.y() << " -> " << ix << " : " << iy << std::endl;
    //    for (int32_t xi = -1; xi <= 1; ++xi) {
    //        for (int32_t yi = -1; yi <= 1; ++yi) {
    //            const auto& cell = getCell(ix + xi, iy + yi);
    //            for (auto j : cell) {
    //                auto& pj = particles[j];
    //                vec r = pj.pos - p.pos;
    //                if (r.squaredNorm() <= 0.99 * 2.0 * 2.0 * packing_2D * packing_2D) {
    //                    emit = false;
    //                    pj.vel = vec(speed2, 0.0);
    //                    pj.angularVelocity = 0.0;
    //                }
    //            }
    //        }
    //    }


    //    //p.vel = vec(speed, 0);
    //    //if (emit) {
    //    //    p.uid = uidCounter++;
    //    //    particles.push_back(p);
    //    //    particlesDFSPH.push_back(dfsphState{});
    //    //}
    //}

    ////    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    ////std::vector<Particle> filteredParticles;
    ////filteredParticles.reserve(particles.size());
    ////for (auto& p : particles) {
    ////}

  scalar vMax = 0.0;
  //static auto& dt = ParameterManager::instance().get<scalar>("sim.dt");
  static auto &damping = ParameterManager::instance().get<scalar>("props.damping");
  for (int32_t i = 0; i < particles.size(); ++i) {
      auto& p = particles[i];
      p.vel += dt * p.accel;
      vMax = std::max(vMax, p.vel.norm());
  }


#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
      auto& p = particles[i];
    //if (p.pos.x() < domainEpsilon * 2.0)
     //   p.vel = vec(speed, 0.0);

    auto pos = p.pos;
    int32_t it = 0;
    for (auto Poly : obstacles) {
        it++;
        //if (it != 3) continue;
        scalar d = sdpoly(Poly, pos);
        scalar dx = sdpoly(Poly, pos + vec(0.001, 0.0));
        scalar dy = sdpoly(Poly, pos + vec(0.0, 0.001));
        vec grad = vec(dx - d, dy - d) / 0.001;
        vec n = grad.normalized();
        vec vortho = (p.vel.dot(n)) * n;
        vec vtangent = p.vel - vortho;
        auto nd = std::max(d, 0.0) / (scale);
        if (d < scale) {
            vtangent *= nd;
            //vtangent *= 0.0;
            //p.vel = vortho + vtangent;
        }
    }

    // shifting starts here

    //int32_t N_i = 0;
    //scalar r_i = 0.0;
    //vec X_i(0., 0.);
    //scalar v_local_max = 0.0, vel_gradient = 0.0;
    //for (int32_t j : p.neighbors) {
    //    auto& pj = particles[j];
    //    vel_gradient += area / pj.rho * (pj.vel - p.vel).dot(gradW(p.pos, pj.pos));
    //    v_local_max = std::max(v_local_max, pj.vel.norm());
    //    if (i == j)
    //        continue;
    //    scalar x_ij = (p.pos - pj.pos).norm();
    //    r_i += x_ij;
    //    N_i++;
    //    X_i += (p.pos - pj.pos) / (x_ij * x_ij * x_ij);
    //}
    //r_i /= N_i;
    //vec X_irr = X_i * r_i * r_i;

    //scalar C = 0.004;
    //scalar a = vMax * dt * C;
    //scalar b = v_local_max * dt * C;
    //scalar d = r_i * 0.001 / dt;
    //scalar f = std::max(b, d);
    //scalar Cdash = std::min(a, f);
    //vec delta_ii = Cdash * X_irr;
    //scalar delta_ii_l = delta_ii.norm();
    //if (N_i < 5)
    //    delta_ii = vec(0.0, 0.0);
    //if(p.pos.x() > 2. * domainEpsilon && p.pos.x() < domainWidth - 2. * domainEpsilon)
    //if (delta_ii_l < 0.75 * support) {
    //    p.pos += delta_ii;
    //    p.vel = p.vel + delta_ii * vel_gradient;
    //}

    scalar k0 = W(vec(0, 0), vec(packing_2D, 0));
    vec X_i(0., 0.);
    scalar v_local_max = 0.0, vel_gradient = 0.0;
    for (int32_t j : p.neighbors) {
        auto& pj = particles[j];
        scalar R = 0.2;
        scalar X_ij = 0.2 * std::pow((W(p.pos, pj.pos) / k0), 4.0);
        vel_gradient += area / pj.rho * (pj.vel - p.vel).dot(gradW(p.pos, pj.pos));
        X_i += area / pj.rho * (1 + X_ij) * gradW(p.pos, pj.pos) * scale * scale;
        v_local_max = std::max(v_local_max, pj.vel.norm());
    }

    scalar C = 250.0; // speed of sound in m/s

    auto Ma = vMax / C;
    auto kCFL = C / scale * dt;

    X_i *= -2.0 * kCFL * Ma;
    auto X_il = X_i.norm();
    auto limit = 0.1 * vMax * dt;

    if (X_il > limit)
        X_i = X_i / X_il * limit;


    // shifting ends here

    //if (delta_ii_l < 0.75f * math::unit_get<4>(p))
    //    arrays.position.second[i] = p + delta_ii;
    //else
    //    arrays.position.second[i] = p;
    //if (delta_ii_l < 0.75f * math::unit_get<4>(p))
    //    arrays.velocity.first[i] = v + delta_ii * vel_gradient;
    //else
    //    arrays.velocity.first[i] = v;

    //std::cout << "[" << p.pos.x() << ", " << p.pos.y() << "] + [" << (dt * p.vel.x()) << ", " << dt * p.vel.y() << "] + [" << X_i.x() << ", " << X_i.y() << "]" << " @ " << kCFL << " : " << Ma << std::endl;


    //auto dist = DBL_MAX;
    //auto x = p.pos;
    //boundaryFunc(p.pos, [&dist](auto bpos, auto d, auto k, auto gk, auto triangle) { dist = std::min(dist,d); });
    //if (p.pos.x() > domainWidth - 3. * domainEpsilon || p.pos.x() < 2. * domainEpsilon)
    //    p.pos += dt * p.vel;
    //else if (dist < scale) {
    //    boundaryFunc(p.pos, [&X_i, &x, &dist](auto bpos, auto d, auto k, auto gk, auto triangle) {
    //        vec n = x - bpos;
    //        if(n.norm() > 1e-5)
    //            n = n / n.norm();

    //        vec vortho = (X_i.dot(n)) * n;
    //        vec vtangent = X_i - vortho;
    //        auto nd = std::max(d, 0.0) / (scale);
    //        if (d < scale) {
    //            vortho *= 0.0;
    //            X_i = vortho + vtangent;
    //        }
    //        });
    //    p.pos += dt * p.vel;
    //    p.vel = p.vel;
    //}
    //else {
    //    p.pos += dt * p.vel + X_i;
    //    p.vel = p.vel + X_i * vel_gradient;
    //}

    p.pos += dt * p.vel;
    p.vel = p.vel * (1. - damping);



    //auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
    //auto t = 12.5 / domainScale;

    //vec begin(domainEpsilon + 10 / domainScale, m - t);
    //vec end(domainWidth - 10 / domainScale, m + t);
    //vec mid = (end + begin) * 0.5;
    //mid.x() -= 25.0 / domainScale;

    //auto thickness = 1.0 / domainScale;
    //auto c = 4. / domainScale;

    //vec center = mid;

    //vec dist = p.pos - center;
    //scalar d = dist.norm() - 3.0 / domainScale;
    //auto grad = dist.normalized();
    //vec n = grad.normalized();
    //vec vortho = (p.vel.dot(n)) * n;
    //vec vtangent = p.vel - vortho;
    //auto nd = std::max(d, 0.0) / (scale);
    //if (d <scale) {
    //    vtangent *= nd;
    //    //vtangent *= 0.0;
    //    p.vel = vortho + vtangent;
    //}


	//if (p.vel.norm() > 1e2)
	//	printParticle(i);
  }
  simulationTime += dt;
  static auto& t = ParameterManager::instance().get<scalar>("sim.time");
  static auto& vMaxr = ParameterManager::instance().get<scalar>("ptcl.maxVelocity");
  vMaxr = vMax;
  t = simulationTime;
  static auto& dtmin = ParameterManager::instance().get<scalar>("sim.minDt");
  static auto& dtmax = ParameterManager::instance().get<scalar>("sim.maxDt");
  dt = std::clamp(0.4 / vMax, dtmin, dtmax);
  static auto& dtr = ParameterManager::instance().get<scalar>("sim.dt");
  dtr = dt;
}
#include <iomanip>

void computeAlpha(bool density = true) {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
    vec kernelSum1(0, 0);
    scalar kernelSum2 = 0.0;
    if (density)
      boundaryFunc(pi.pos, [&pi, &kernelSum1, &kernelSum2](auto bpos, auto d, auto k, auto gk, auto triangle) { kernelSum1 += gk; });
    for (int32_t j : pi.neighbors) {
      auto &pj = particles[j];
      auto &dj = particlesDFSPH[j];
      vec kernel = gradW(pi.pos, pj.pos);
      kernelSum1 += dj.area * kernel;
      kernelSum2 += dj.area * dj.area / mass * kernel.dot(kernel);
    }
    di.alpha = -dt * dt * di.area / mass * kernelSum1.dot(kernelSum1) - dt * dt * di.area * kernelSum2;
	if (std::abs(di.alpha) > 1.0) {
		boundaryFunc(pi.pos, [&pi, &kernelSum1, &kernelSum2](auto pb, auto d, auto k, auto gk, auto triangle) { 
			kernelSum1 += gk; 
			if (triangle)
				std::cout << "Triangle:" << std::endl;
			std::cout << std::setprecision(12) << pi.pos.x() << " : " << std::setprecision(12) << pi.pos.y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << d << ", " << k << " ==> " << gk.x() << " : " << gk.y() << std::endl;
			});

	}
  }
}
void computeSourceTerm(bool density = true) {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
    scalar sourceTerm = density ? scalar(1.0) - pi.rho : scalar(0.0);
    if (density)
      boundaryFunc(pi.pos, [&di, &sourceTerm](auto bpos, auto d, auto k, auto gk, auto triangle) {
        sourceTerm = sourceTerm - dt * di.vel.dot(gk);
      });
    for (int32_t j : pi.neighbors) {
      auto &pj = particles[j];
      auto &dj = particlesDFSPH[j];
      sourceTerm -= dt * dj.area * (di.vel - dj.vel).dot(gradW(pi.pos, pj.pos));
    }
    di.source = sourceTerm;
    di.pressure2 = (scalar)0.0;
  }
}
void computeAcceleration(bool density = true) {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
    vec kernelSum((scalar)0.0, (scalar)0.0);
    for (int32_t j : pi.neighbors) {
      auto &pj = particles[j];
      auto &dj = particlesDFSPH[j];
      kernelSum += -mass * (di.pressure2 / power(pi.rho * rho0, 2) + dj.pressure2 / power(pj.rho * rho0, 2)) *
                   gradW(pi.pos, pj.pos);
    }
    di.accel += kernelSum;
    di.pressure1 = di.pressure2;
  }
}
void updatePressure(bool density = true) {
    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
    scalar kernelSum = (scalar)0.0;
    if (density)
      boundaryFunc(pi.pos,
                   [&di, &kernelSum](auto bpos, auto d, auto k, auto gk, auto triangle) { kernelSum += dt * dt * di.accel.dot(gk); });
    for (int32_t j : pi.neighbors) {
      auto &pj = particles[j];
      auto &dj = particlesDFSPH[j];
      kernelSum += dt * dt * dj.area * (di.accel - dj.accel).dot(gradW(pi.pos, pj.pos));
    }
    scalar omega = (scalar)0.5;
    scalar pressure = di.pressure1 + omega / di.alpha * (di.source - kernelSum);
    //if(di.pressureBoundary == 0.0)
    if(backgroundPressureSwitch)
        pressure += 1.0 * speed * speed * rho0;
    pressure = density ? std::max(pressure, (scalar)0.0) : pressure;
    scalar residual = kernelSum - di.source;
    if (::abs(di.alpha) < 1e-25 || pressure != pressure || pressure > 1e25)
      pressure = residual = 0.0;
    //pressure = std::max(pressure, speed * speed * rho0);
    di.pressure2 = pressure;
    di.dpdt = std::max(residual, (scalar)-0.001) * area;
    di.rhoStar = std::max(residual, (scalar)-0.001) * area;
  }
}
scalar calculateBoundaryPressureMLS(int32_t i, vec pb, bool density = true) {
  vec vecSum(0, 0), d_bar(0, 0);
  matrix M = matrix::Zero();
  scalar sumA = (scalar)0.0, sumB = (scalar)0.0, d_sum = (scalar)0.0;

  auto [ix, iy] = getCellIdx(pb.x(), pb.y());
  for (int32_t xi = -1; xi <= 1; ++xi) {
    for (int32_t yi = -1; yi <= 1; ++yi) {
      const auto &cell = getCell(ix +xi, iy +yi);
      for (auto j : cell) {
        auto &pj = particles[j];
        vec r = pj.pos - pb;
        if (r.squaredNorm() <= support * support) {
          scalar fac = area * W(pj.pos, pb);
          d_bar += pj.pos * fac;
          d_sum += fac;
        }
      }
    }
  }
  d_bar /= d_sum;
  vec x_b = pb - d_bar;

  for (int32_t xi = -1; xi <= 1; ++xi) {
    for (int32_t yi = -1; yi <= 1; ++yi) {
        const auto& cell = getCell(ix + xi, iy + yi);
      for (auto j : cell) {
        auto &pj = particles[j];
        auto &dj = particlesDFSPH[j];
        vec r = pj.pos - pb;
        if (r.squaredNorm() <= (scalar)1.0) {
          scalar Wbbf = W(pj.pos, pb);
          vec pjb = pj.pos - d_bar;
          M += (pjb) * (pjb).transpose();
          vecSum += pjb * dj.pressure2 * area * Wbbf;
          sumA += dj.pressure2 * area * Wbbf;
          sumB += area * Wbbf;
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
void computeBoundaryPressure(bool density = true) {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
	di.pressureBoundary = 0.0;
    di.accel = vec(0, 0);
    if (!density)
      continue;
    boundaryFunc(pi.pos, [&](auto pb, auto d, auto k, auto gk, auto triangle) {
      scalar pressure = calculateBoundaryPressureMLS(i, pb, density);
      scalar fluidPressure = density ? std::max((scalar)0.0, di.pressure2) : di.pressure2;
      scalar boundaryPressure = density ? std::max((scalar)0.0, pressure) : pressure;
	 // boundaryPressure = fluidPressure;
	  di.pressureBoundary += boundaryPressure;
      di.accel +=
          -scalar(1.0) * rho0 * (fluidPressure / power(pi.rho * rho0, 2) + boundaryPressure / (rho0 * rho0)) * gk;
	  if (boundaryPressure > 1e12 && triangle) {
		  std::cout << pi.pos.x() << " : " << pi.pos.y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << d << ", " << k << " ==> " << gk.x() << " : " << gk.y() << std::endl;
		  printParticle(i);
	  }
    });
  }
}
void predictVelocity(bool density = true) {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
    di.vel = pi.vel + dt * pi.accel;
    di.area = area / pi.rho;
  }
}
void updateVelocity(bool density = true) {
#pragma omp parallel for
  for (int32_t i = 0; i < particles.size(); ++i) {
    auto &pi = particles[i];
    auto &di = particlesDFSPH[i];
    pi.accel += di.accel;
    di.vel += dt * di.accel;
  }
  
}
int32_t divergenceSolve() {
	//return 0;
  predictVelocity();
  computeAlpha(false);
  computeSourceTerm(false);
  //scalar error = (scalar)0.0;
  //int32_t counter = 0;
  scalar totalArea = area * (scalar)particles.size();
  static auto& limit = ParameterManager::instance().get<scalar>("dfsph.divergenceEta");
  static auto& error = ParameterManager::instance().get<scalar>("dfsph.divergenceError");
  static auto& counter = ParameterManager::instance().get<int32_t>("dfsph.divergenceIterations");
  counter = 0;
  do {
    computeBoundaryPressure(false);
    computeAcceleration(false);
    updatePressure(false);
    error = (scalar)0.0;
    for (auto di : particlesDFSPH)
      error += di.dpdt;
    error /= totalArea;
    //std::cout << "Divergence: " << counter << " -> " << error << std::endl;
  } while (counter++ < 2 || (error > (scalar)limit && counter < 256));
  computeBoundaryPressure(false);
  computeAcceleration(false);
  updateVelocity(false);
  return counter;
}
int32_t densitySolve() {
  predictVelocity();
  computeAlpha(true);
  computeSourceTerm(true);
  //scalar error = (scalar)0.0;
  //int32_t counter = 0;
  for (auto& dp : particlesDFSPH) {
      dp.pressure1 = dp.pressure2 = 0.5 * dp.pressure1;
  }
  scalar totalArea = area * (scalar)particles.size();
  static auto& limit = ParameterManager::instance().get<scalar>("dfsph.densityEta");
  static auto& error = ParameterManager::instance().get<scalar>("dfsph.densityError");
  static auto& counter = ParameterManager::instance().get<int32_t>("dfsph.densityIterations");
  counter = 0;
  do {
    computeBoundaryPressure(true);
    computeAcceleration(true);
    updatePressure(true);
    error = (scalar)0.0;
    for (auto di : particlesDFSPH)
      error += di.rhoStar;
    error /= totalArea;
   // std::cout << "Density: " << counter << " -> " << error << std::endl;
  } while (counter++ < 3 || (error > (scalar)limit && counter < 32));
  computeBoundaryPressure(true);
  computeAcceleration(true);
  updateVelocity();


  for (auto& dp : particlesDFSPH) {
      dp.pressurePrevious = dp.pressure1;
  }

  return counter;
}

void resetFrame() {
    
  for (auto &v : cellArray)
    v.clear();
  for (auto &p : particles)
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

  if (particlesDFSPH.size() != particles.size())
    particlesDFSPH.resize(particles.size());
  for (auto &d : particlesDFSPH)
    d.reset();
  static auto& numPtcls = ParameterManager::instance().get<std::size_t>("props.numptcls");
  numPtcls = particles.size();
} 

void emitParticles() {
    if(!inletSwitch)
    return;
    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
    //auto dt = ParameterManager::instance().get<scalar>("sim.dt");
    static scalar offset = -dt * speed;
    offset = std::fmod(offset + dt * speed, 2.0 * packing_2D);

    auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
    auto t = 12.5 / domainScale;
    auto eps = scale/ sqrt(2) * 0.2;
    auto particlesGen = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon+ spacing_2D), vec(domainEpsilon + 5.0/domainScale, domainHeight- domainEpsilon - spacing_2D));
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
                    if (r.squaredNorm() <= 0.5 *2.0 * 2.0 * packing_2D * packing_2D) {
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
#include <random>
#include <iterator>

void initializeSPH(int32_t scene) {
    ParameterManager::instance().newParameter("ray.origin", vec(50,25), { .constant = false});
    ParameterManager::instance().newParameter("ray.target", vec(5,5), { .constant = false });
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
    ParameterManager::instance().newParameter("colorMap.buffer", std::string("angularVelocity"), { .constant = false ,
        .presets = std::vector<detail::iAny>{
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
            std::string("pressureAcceleration")}
        });
    ParameterManager::instance().newParameter("colorMap.vectorMode", std::string("magnitude"), { .constant = false ,
        .presets = std::vector<detail::iAny>{
            std::string("magnitude"),
            std::string("x"),
            std::string("y")}
        });
    ParameterManager::instance().newParameter("colorMap.map", 0, { .constant = false , .range = Range{0,3} });
    ParameterManager::instance().newParameter("colorMap.limit", true, { .constant = false });
    ParameterManager::instance().newParameter("colorMap.auto", true, { .constant = false});
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
    ParameterManager::instance().newParameter("ptcl.viscosityConstant", 0.01, { .constant = false, .range = Range{0.01, 0.05} });

    ParameterManager::instance().newParameter("dfsph.densityEta", scalar(0.001), { .constant = true });
    ParameterManager::instance().newParameter("dfsph.divergenceEta", scalar(0.001), { .constant = true });
    ParameterManager::instance().newParameter("dfsph.densityIterations", 0, { .constant = true });
    ParameterManager::instance().newParameter("dfsph.divergenceIterations", 0, { .constant = true });
    ParameterManager::instance().newParameter("dfsph.densityError", scalar(0), { .constant = true });
    ParameterManager::instance().newParameter("dfsph.divergenceError", scalar(0), { .constant = true });


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
    vec P1(domainEpsilon, domainEpsilon);
    vec P2(domainWidth - domainEpsilon, domainEpsilon);
    vec P3(domainWidth - domainEpsilon, domainHeight - domainEpsilon);
    vec P4(domainEpsilon, domainHeight - domainEpsilon);

    polygon.push_back(P2);
    polygon.push_back(P3);
    polygon.push_back(P4);
    polygon.push_back(P1);
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

    vec begin(domainEpsilon+10 / domainScale, m-t);
    vec end(domainWidth -20 / domainScale, m+t);
    vec mid = (end + begin) * 0.5;
    mid.x() -= 25.0 / domainScale;
    mid.y() -= 0.1 / domainScale;

    auto thickness = 1.0 / domainScale;
    auto c = 3.5 / domainScale;


    std::vector<vec> Upper{{begin.x(), begin.y()},
                           {end.x(), begin.y()},
                           {end.x(), begin.y() - thickness},
                           {begin.x(), begin.y() - thickness}};
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
    std::cout << "New Gap: " << left << " [ " << ::ceil(left) * packing_2D << " + " << domainEpsilon << " + " << 2 * spacing_2D << " ]"<< std::endl;
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


    //obstacles.push_back(Upper);
    //obstacles.push_back(Lower);
    //obstacles.push_back(Obs);
    switch (bc) {
    case boundaryConfig::Box:
        obstacles.push_back(ObsBox); break;
    case boundaryConfig::CenterBox:
        obstacles.push_back(ObsBoxCenter); break;
    case boundaryConfig::ObstacleBT:
        obstacles.push_back(ObsB);
        obstacles.push_back(ObsT); break;
    case boundaryConfig::ObstacleBTLow:
        obstacles.push_back(ObsB2);
        obstacles.push_back(ObsT2); break;
    case boundaryConfig::Trapezoid:
        obstacles.push_back(Obs); break;
    }
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
        switch (pc) {
        case particleConfig::Domain:particles5 = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainWidth - domainEpsilon - spacing_2D, domainHeight - domainEpsilon - spacing_2D)); break;
        case particleConfig::DamBreak:particles5 = genParticles(vec(domainEpsilon + spacing_2D, domainEpsilon + spacing_2D), vec(domainEpsilon + 20.0 / domainScale, domainEpsilon + 37.5 / domainScale)); break;
        }

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
