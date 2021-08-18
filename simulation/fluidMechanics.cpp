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

void density() {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            pi.rho += area * W(pi.pos, pj.pos);
        }
        //pi.rho = std::max(pi.rho, 0.5);
        boundaryFunc(pi, [&pi, i](auto bpos, auto d, auto k, auto gk, auto triangle) {
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
    std::vector<vec> dvdt(particles.size(), vec(0.0, 0.0));
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        scalar voriticity = 0.0;
        vec vterm(0.0, 0.0);

        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto grad = gradW(pi.pos, pj.pos);
            auto vel = pi.vel - pj.vel;
            auto vor = pi.angularVelocity - pj.angularVelocity;
            auto term = area / pj.rho * (vel.x() * grad.y() - vel.y() * grad.x());

            voriticity += term;
            vterm.x() += area / pj.rho * (-vor * grad.y());
            vterm.y() += area / pj.rho * (vor * grad.x());


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
    if (gravitySwitch)
        for (auto& p : particles)
            p.accel += gravity;
}

void XSPH() {
    static auto& viscosityConstant = ParameterManager::instance().get<scalar>("ptcl.viscosityConstant");
    std::vector<vec> tempV;
    //#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        tempV.push_back(pi.vel);
        for (auto& j : pi.neighbors) {
            auto& pj = particles[j];
            tempV[i] += viscosityConstant * (pi.rho + pj.rho) * area / (pi.rho + pj.rho) * scalar(2.0) * W(pi.pos, pj.pos) * (pj.vel - pi.vel);
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
        if (outletSwitch)
            if (p.pos.x() > domainWidth - 2. * domainEpsilon) {
                p.vel = vec(speed2, 0.0);
                p.angularVelocity = 0.0;
                p.accel = vec(0.0, 0.0);
            }
        if (inletSwitch)
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
    static auto& damping = ParameterManager::instance().get<scalar>("props.damping");
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
                //vtangent *= nd;
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
    dt = std::clamp(0.4 * support / vMax, dtmin, dtmax);
    static auto& dtr = ParameterManager::instance().get<scalar>("sim.dt");
    dtr = dt;
}
#include <iomanip>

void computeAlpha(bool density) {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        vec kernelSum1(0, 0);
        scalar kernelSum2 = 0.0;
        if (density)
            boundaryFunc(pi, [&pi, &kernelSum1, &kernelSum2](auto bpos, auto d, auto k, auto gk, auto triangle) { kernelSum1 += gk; });
        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto& dj = particlesDFSPH[j];
            vec kernel = gradW(pi.pos, pj.pos);
            kernelSum1 += dj.area * kernel;
            kernelSum2 += dj.area * dj.area / mass * kernel.dot(kernel);
        }
        di.alpha = -dt * dt * di.area / mass * kernelSum1.dot(kernelSum1) - dt * dt * di.area * kernelSum2;
        if (std::abs(di.alpha) > 1.0) {
            boundaryFunc(pi, [&pi, &kernelSum1, &kernelSum2](auto pb, auto d, auto k, auto gk, auto triangle) {
                kernelSum1 += gk;
                if (triangle)
                    std::cout << "Triangle:" << std::endl;
                std::cout << std::setprecision(12) << pi.pos.x() << " : " << std::setprecision(12) << pi.pos.y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << d << ", " << k << " ==> " << gk.x() << " : " << gk.y() << std::endl;
                });

        }
    }
}
void computeSourceTerm(bool density) {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        scalar sourceTerm = density ? scalar(1.0) - pi.rho : scalar(0.0);
        if (density)
            boundaryFunc(pi, [&di, &sourceTerm](auto bpos, auto d, auto k, auto gk, auto triangle) {
            sourceTerm = sourceTerm - dt * di.vel.dot(gk);
                });
        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto& dj = particlesDFSPH[j];
            sourceTerm -= dt * dj.area * (di.vel - dj.vel).dot(gradW(pi.pos, pj.pos));
        }
        di.source = sourceTerm;
        di.pressure2 = (scalar)0.0;
    }
}
void computeAcceleration(bool density) {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        vec kernelSum((scalar)0.0, (scalar)0.0);
        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto& dj = particlesDFSPH[j];
            kernelSum += -mass * (di.pressure2 / power(pi.rho * rho0, 2) + dj.pressure2 / power(pj.rho * rho0, 2)) *
                gradW(pi.pos, pj.pos);
        }
        di.accel += kernelSum;
        di.pressure1 = di.pressure2;
    }
}
void updatePressure(bool density) {
    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        scalar kernelSum = (scalar)0.0;
        if (density)
            boundaryFunc(pi,
                [&di, &kernelSum](auto bpos, auto d, auto k, auto gk, auto triangle) { kernelSum += dt * dt * di.accel.dot(gk); });
        for (int32_t j : pi.neighbors) {
            auto& pj = particles[j];
            auto& dj = particlesDFSPH[j];
            kernelSum += dt * dt * dj.area * (di.accel - dj.accel).dot(gradW(pi.pos, pj.pos));
        }
        scalar omega = (scalar)0.5;
        scalar pressure = di.pressure1 + omega / di.alpha * (di.source - kernelSum);
        //if(di.pressureBoundary == 0.0)
        if (backgroundPressureSwitch)
            pressure += 1.0 * speed * speed * rho0;
        pressure = density ? std::max(pressure, (scalar)0.0) : pressure;
        scalar residual = kernelSum - di.source;
        if (::abs(di.alpha) < 1e-25 || pressure != pressure || pressure > 1e25)
            pressure = residual = 0.0;
        //pressure = std::max(pressure, speed * speed * rho0);
        di.pressure2 = pressure;
        di.dpdt = std::max(residual, (scalar)-0.001) * area;
        di.rhoStar = residual * area;
    }
}
scalar calculateBoundaryPressureMLS(int32_t i, vec pb, bool density) {
    vec vecSum(0, 0), d_bar(0, 0);
    matrix M = matrix::Zero();
    scalar sumA = (scalar)0.0, sumB = (scalar)0.0, d_sum = (scalar)0.0;

    int32_t ii = 0;
    auto [ix, iy] = getCellIdx(pb.x(), pb.y());
    for (int32_t xi = -1; xi <= 1; ++xi) {
        for (int32_t yi = -1; yi <= 1; ++yi) {
            if (xi + ix < 0 || xi + ix >= cellsX) continue;
            if (yi + iy < 0 || yi + iy >= cellsY) continue;
            const auto& cell = getCell(ix + xi, iy + yi);
            for (auto j : cell) {
                auto& pj = particles[j];
                vec r = pj.pos - pb;
                if (r.squaredNorm() <= support * support) {
                    ii++;
                    scalar fac = area * W(pj.pos, pb);
                    d_bar += pj.pos * fac;
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
                auto& pj = particles[j];
                auto& dj = particlesDFSPH[j];
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
    std::cout << std::defaultfloat;
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
void computeBoundaryTrianglePressure(bool density) {
    static auto& baryPressure = ParameterManager::instance().get<bool>("sim.barycentricPressure");
    if (!baryPressure) return;
#pragma omp parallel for
    for (int32_t i = 0; i < triangles.size(); ++i) {
        const auto [t0, t1, t2] = triangles[i];
        auto f0 = calculateBoundaryPressureMLS(-1, t0, true);
        auto f1 = calculateBoundaryPressureMLS(-1, t1, true);
        auto f2 = calculateBoundaryPressureMLS(-1, t2, true);
        trianglePressures[i] = std::make_tuple(f0, f1, f2);
    }

}
void computeBoundaryPressure(bool density) {
    static auto& baryPressure = ParameterManager::instance().get<bool>("sim.barycentricPressure");
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        di.pressureBoundary = 0.0;
        di.accel = vec(0, 0);
        if (!density)
            continue;

        if (baryPressure) {
            for (auto ti : pi.neighborTriangles) {
                const auto& tri = triangles[ti];
                auto [f0, f1, f2] = trianglePressures[ti];
                scalar fluidPressure = density ? std::max((scalar)0.0, di.pressure2) : di.pressure2;
                auto [pb, d] = closestPointTriangle(pi.pos, tri);
                scalar pressure = calculateBoundaryPressureMLS(0, pb, density);
                //f0 = f1 = f2 = pressure;
                //f0 = std::max(f0,std::max(f1, f2));
                //f1 = f2 = f0;
                auto [hit, pb2, d2, k, gk] = interactTriangleBaryCentric(pi.pos, tri, pi.rho * rho0, fluidPressure, f0, f1, f2);
                //auto [hit2, pb3, d3, k3, gk3] = interactTriangle(pi.pos, tri);


                if (hit) {
                    //auto [hit, pb, d, k2, gk2] = interactTriangleBaryCentric(pi.pos, tri, pi.rho * rho0, fluidPressure, pressure, pressure, pressure);
                    //std::cout << "Hit Triangle!!!!\n\n\n\n";
                    //if (std::abs(gk.norm() - gk2.norm())>1e-2)
                    //{
                    //    std::cout << "\n";
                    //    std::cout << std::setprecision(5) << f0 << " " << f1 << " " << f2 << " -> " << pressure << " => " << k << " @ [ " << gk.x() << " " << gk.y() << "] : " << k3 << " @ [ " << gk3.x() << " " << gk3.y() << "] : ";
                    //    std::cout << k2 << " @ [ " << gk2.x() << " " << gk2.y() << "] \n";
                    //}
                    di.accel += -scalar(1.0) * gk;
                }
                //if (f0 > 1e12 || f1 > 1e12 || f2 > 1e12) {
                //    std::cout << std::defaultfloat;
                //    std::cout << "\n\n";
                //    std::cout << f0 << " " << f1 << " " << f2 << " - " << fluidPressure << std::endl;
                //    std::cout << pi.pos.x() << " : " << pi.pos.y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << d << ", " << k << " ==> " << gk.x() << " : " << gk.y() << std::endl;
                //    printParticle(i);
                //}

                di.pressureBoundary += f0 + f1 + f2;
            }
        }
        else {
            boundaryFunc(pi, [&](auto pb, auto d, auto k, auto gk, auto triangle) {
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
}
void predictVelocity(bool density) {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        di.vel = pi.vel + dt * pi.accel;
        di.area = area / pi.rho;
    }
}
void updateVelocity(bool density) {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        pi.accel += di.accel;
        di.vel += dt * di.accel;
    }
}
void predictVelocityPCI(bool density) {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        di.vel = pi.vel + dt * pi.accel + dt * di.accel;
        di.pos = pi.pos + dt * di.vel;
        di.area = area / pi.rho;
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
        //computeBoundaryTrianglePressure(false);
        computeBoundaryPressure(false);
        computeAcceleration(false);
        updatePressure(false);
        error = (scalar)0.0;
        for (auto di : particlesDFSPH)
            error += di.dpdt;
        error /= totalArea;
        //std::cout << "Divergence: " << counter << " -> " << error << std::endl;
    } while (counter++ < 3 || (error > (scalar)limit && counter < 3));
    computeBoundaryPressure(false);
    computeAcceleration(false);
    updateVelocity(false);
    return counter;
}
void predictDensity() {
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];
        auto [ix, iy] = getCellIdx(di.pos.x(), di.pos.y());
        di.rhoStar = 0.;
        for (int32_t xi = -2; xi <= 2; ++xi) {
            for (int32_t yi = -2; yi <= 2; ++yi) {
                const auto& cell = getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    auto& pj = particlesDFSPH[j];
                    vec r = pj.pos - di.pos;
                    if (r.squaredNorm() <= support * support)
                        di.rhoStar += area * W(pi.pos, pj.pos);
                }
            }
        }
        //pi.rho = std::max(pi.rho, 0.5);
        boundaryFunc(di.pos, [&di, i](auto bpos, auto d, auto k, auto gk, auto triangle) {
            //std::cout << i << " - " << pi.pos.x() << " : " << pi.pos.y() << " -> " << bpos.x() << " : " << bpos.y() << ", " << d << ", " << k << ", " << gk.x() << " : " << gk.y() << std::endl;
            di.rhoStar += k; });
        di.dpdt = di.rhoStar - 1.;
    }
}
void updatePressurePCI() {
    auto speed = ParameterManager::instance().get<scalar>("sim.inletSpeed");
#pragma omp parallel for
    for (int32_t i = 0; i < particles.size(); ++i) {
        auto& pi = particles[i];
        auto& di = particlesDFSPH[i];

        auto delta = 0x1.41719825bf00dp-32;
        auto beta = dt * dt * area * area * 2.;
        auto pressure = rho0 * (di.rhoStar - 1.0) * delta / beta;
        pressure = std::max(0., pressure);
        di.pressure1 += pressure;
        di.pressure2 = di.pressure1;
    }

}
int32_t PCISPH() {
    static auto& limit = ParameterManager::instance().get<scalar>("dfsph.densityEta");
    static auto& error = ParameterManager::instance().get<scalar>("dfsph.densityError");
    static auto& counter = ParameterManager::instance().get<int32_t>("dfsph.densityIterations");
    scalar totalArea = (scalar)particles.size();
    counter = 0;
    for (auto& dp : particlesDFSPH) {
        dp.pressure1 = dp.pressure2 = 0.;
    }

    do {
        computeBoundaryPressure(true);
        computeAcceleration(true);
        predictVelocityPCI(true);
        predictDensity();

        updatePressurePCI();




        error = (scalar)0.0;
        for (auto di : particlesDFSPH)
            error += std::max(-0.005, di.dpdt);
        error /= totalArea;
        //std::cout << "Density: " << counter << " -> " << std::defaultfloat << error << std::endl;
    } while (counter++ < 3 || (error > (scalar)limit * 10. && counter < 32));
    computeBoundaryPressure(true);
    computeAcceleration(true);
    updateVelocity();

    return counter;
}

int32_t densitySolve() {
    // return PCISPH();
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
        computeBoundaryTrianglePressure(true);
        computeBoundaryPressure(true);
        computeAcceleration(true);
        updatePressure(true);
        error = (scalar)0.0;
        for (auto di : particlesDFSPH)
            error += std::max(-0.001 * area, di.rhoStar);
        error /= totalArea;
        // std::cout << "Density: " << counter << " -> " << error << std::endl;
    } while (counter++ < 3 || (error > (scalar)limit && counter < 256));
    computeBoundaryPressure(true);
    computeAcceleration(true);
    updateVelocity();


    for (auto& dp : particlesDFSPH) {
        dp.pressurePrevious = dp.pressure1;
    }

    return counter;
}