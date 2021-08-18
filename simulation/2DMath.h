#pragma once
#define _USE_MATH_DEFINES
#include "Math.h"
#include "SPH.h"
#include <math.h>
// this set of power functions is used to avoid any potential inaccuracies of calling a normal pow function with integer
// powers and always uses a simple iterative multiplication of the base.
inline auto power(std::complex<scalar> v, int32_t p) {
    if (p == 0)
        return std::complex<scalar>(scalar(0.0), scalar(0.0));
    std::complex<scalar> s = v;
    for (int32_t i = 1; i < p; ++i)
        s *= v;
    return s;
}
inline auto power(scalar v, int32_t p) {
    if (p == 0)
        return scalar(0.0);
    scalar s = v;
    for (int32_t i = 1; i < p; ++i)
        s *= v;
    return s;
}
// for some functions it is easier to describe complex values as 2.0 + 1.0_i instead of requiring std::complex(2,1)
inline constexpr std::complex<scalar> operator"" _i(long double d) {
    return std::complex<scalar>(scalar(0.0), static_cast<scalar>(d));
}
inline constexpr std::complex<scalar> operator"" _i(unsigned long long d) {
    return std::complex<scalar>(scalar(0.0), static_cast<scalar>(d));
}
// checks if a given scalar is zero, i.e. abs(scalar) < epsilon
inline auto isZero(scalar v) { return std::abs(v) < epsilon; }
// returns the sign of the given scalar. 1 if v > 0, -1 if v < 1.0, 0 else
inline scalar sgn(scalar v) { return v < scalar(-epsilon) ? scalar(-1.0) : (v > scalar(epsilon) ? scalar(1.0) : scalar(0.0)); }
// computes the SVD of a given 2x2 matrix using a closed form solution
inline std::tuple<matrix, matrix, matrix> svd2x2(matrix M) {
    scalar a = M(0, 0);
    scalar a2 = a * a;
    scalar b = M(0, 1);
    scalar b2 = b * b;
    scalar c = M(1, 0);
    scalar c2 = c * c;
    scalar d = M(1, 1);
    scalar d2 = d * d;

    scalar theta = scalar(0.5) * std::atan2(scalar(2.0) * a * c + scalar(2.0) * b * d, a2 + b2 - c2 - d2);
    scalar ct = std::cos(theta);
    scalar st = std::sin(theta);
    matrix U;
    U << ct, -st, st, ct;

    scalar S1 = a2 + b2 + c2 + d2;
    scalar S2 = std::sqrt(power(a2 + b2 - c2 - d2, 2) + scalar(4.0) * power(a * c + b * d, 2));
    scalar sigma1 = std::sqrt((S1 + S2) * scalar(0.5));
    scalar sigma2 = std::sqrt((S1 - S2) * scalar(0.5));
    matrix S;
    S << sigma1, scalar(0.0), scalar(0.0), sigma2;

    scalar phi = scalar(0.5) * std::atan2(scalar(2.0) * a * b + scalar(2.0) * c * d, a2 - b2 + c2 - d2);
    auto cp = std::cos(phi);
    auto sp = std::sin(phi);
    scalar s11 = (a * ct + c * st) * cp + (b * ct + d * st) * sp;
    scalar s22 = (a * st - c * ct) * sp + (-b * st + d * ct) * cp;
    matrix V;
    V << sgn(s11) * cp, -sgn(s22) * sp, sgn(s11)* sp, sgn(s22)* cp;

    return std::make_tuple(U, S, V);
}
// simple cubic spline kernel function. Assumes H = 1.0 for all particles.
inline scalar W(vec a, vec b) {
    vec r = b - a;
    scalar q = r.norm() / support;
    //spline
    //if (q < scalar(0.5))
    //  return kernelConstant * (power(scalar(1.0) - q, 3) - scalar(4.0) * power(scalar(0.5) - q, 3));
    //else if (q <= scalar(1.0))
    //  return kernelConstant * power(scalar(1.0) - q, 3);
    if (q < 1.0)
        return 7.0 / std::numbers::pi / (support * support) * power(1.0 - q, 4) * (1.0 + 4.0 * q);
    // wend 4
    //if (q < 1.0)
    //    return 9.0 / pi / (support * support) * power(1.0 - q, 6) * (1.0 + 6.0 * q + 35.0 / 3.0 * q * q);
    return scalar(0.0);
}
// simple spiky kernel gradient function. Assumes H = 1.0
inline vec gradW(vec a, vec b) {
    vec r = a - b;
    auto rn = r.norm();
    auto q = rn / support;
    if (q < epsilon || q > scalar(1.0))
        return vec(scalar(0.0), scalar(0.0));

    //spline
    //if (q < scalar(0.5))
    //    return -r/rn * gradConstant * (3.0 * power(scalar(1.0) - q, 2) - scalar(12.0) * power(scalar(0.5) - q, 2));
    //else if (q <= scalar(1.0))
    //    return -r / rn * gradConstant *3.0 *  power(scalar(1.0) - q, 2);


    // spiky
    //return -r/rn * scalar(30.0) / scalar(pi * support * support * support) * power(scalar(1.0) - q, 2);

    // wend2
    return -r / rn * 7.0 / std::numbers::pi / (support * support * support) * (20.0 * q * power(1 - q, 3));

    // wend4
    //return -r / rn * 9.0 / pi / (support * support * support) * 56.0 / 3.0 * q * (5.0 * q + 1.0) * power(1 - q, 5);

    // these lines represent a normal cubic spline gradient function, also assuming H = 1.0
    // if (q <= scalar(0.5))
    //  return r * kernelConstant * (-scalar(3.0) * power(scalar(1.0) - q, 2) + scalar(12.0) * power(scalar(0.5) - q, 2));
    // return -scalar(3.0) * r * kernelConstant * power(scalar(1.0) - q, 2);
}
// if POLYFIT is defined a 16th degree polynomial is used to approximate the boundary integral
// instead of using an exact solution.
//#define POLYFIT
// returns the boundary integral for a flat infinite boundary at distance dr
inline scalar boundaryKernelAnalytic(scalar dr) {
    std::complex<scalar> d(-dr, 0.0);
    auto a1 = isZero(dr)
        ? 0.0
        : (+12.0 * std::pow(d, 5) + 80.0 * std::pow(d, 3)) * std::log(std::sqrt(1.0 - std::pow(d, 2)) + 1.0);
    auto a2 = isZero(dr)
        ? 0.0
        : (-12.0 * std::pow(d, 5) - 80.0 * std::pow(d, 3)) * std::log(1.0 - std::sqrt(1.0 - std::pow(d, 2)));
    auto a3 = isZero(dr) ? 0.0
        : (-12.0 * std::pow(d, 5) - 80.0 * std::pow(d, 3)) *
        std::log(std::sqrt(1.0 - 4.0 * std::pow(d, 2)) + 1.0);
    auto a4 = isZero(dr) ? 0.0
        : (+12.0 * std::pow(d, 5) + 80.0 * std::pow(d, 3)) *
        std::log(1.0 - std::sqrt(1.0 - 4.0 * std::pow(d, 2)));
    auto a5 = -13.0 * std::acos(2.0 * d);
    auto a6 = +16.0 * std::acos(d);
    auto a7 = std::sqrt(1.0 - 4.0 * std::pow(d, 2)) * (74.0 * std::pow(d, 3) + 49.0 * d);
    auto a8 = std::sqrt(1.0 - std::pow(d, 2)) * (-136.0 * std::pow(d, 3) - 64.0 * d);
    std::complex<scalar> a = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8;

    auto b1 = -36.0 * std::pow(d, 5) * std::log(std::sqrt(1.0 - 4.0 * std::pow(d, 2)) + 1.0);
    auto b2 = isZero(dr) ? 0.0 : 36.0 * std::pow(d, 5) * std::log(1.0 - std::sqrt(1.0 - 4.0 * power(d, 2)));
    auto b3 = 11.0 * std::acos(2.0 * d);
    auto b4 = -36.0 * std::log(-1.0 + 0.0_i) * std::pow(d, 5);
    auto b5 = -160_i * std::pow(d, 4);
    auto b6 = std::sqrt(1.0 - 4.0 * std::pow(d, 2)) * (62.0 * std::pow(d, 3) - 33.0 * d);
    auto b7 = 80_i * std::pow(d, 2);
    std::complex<scalar> b = b1 + b2 + b3 + b4 + b5 + b6 + b7;
    auto result = ((a + b) / (14.0 * M_PI));
    return result.real();
}
inline scalar boundaryKernelPolyFit(scalar dr) {
    scalar poly = (scalar)1.72460000e-01;
    poly = poly * dr - (scalar)1.13945615e-01;
    poly = poly * dr - (scalar)6.79631508e-01;
    poly = poly * dr + (scalar)4.71685815e-01;
    poly = poly * dr + (scalar)1.07549979e+00;
    poly = poly * dr - (scalar)1.16258333e+00;
    poly = poly * dr - (scalar)8.74915460e-01;
    poly = poly * dr + (scalar)2.30265078e+00;
    poly = poly * dr + (scalar)3.89166005e-01;
    poly = poly * dr - (scalar)3.40417812e+00;
    poly = poly * dr - (scalar)9.29508384e-02;
    poly = poly * dr + (scalar)3.48600557e+00;
    poly = poly * dr + (scalar)1.08432767e-02;
    poly = poly * dr - (scalar)2.44291025e+00;
    poly = poly * dr - (scalar)4.94944643e-04;
    poly = poly * dr + (scalar)1.36322382e+00;
    poly = poly * dr + (scalar)5.00006348e-01;
    return poly;
}
#include <cheb/cheb.h>
inline scalar boundaryKernelChebFit(scalar dr) {
    static cheb::Function chebFn([](cheb::scalar s) {
        return boundaryKernelAnalytic(s);
        });
    return chebFn(dr);
}


#define PRINTVAR(x) {std::cout << #x << "\t=\t"<< x << std::endl;}
inline scalar boundaryKernel(scalar dr) {
    dr = std::clamp(dr, (scalar)-1.0, (scalar)1.0);
    //return boundaryKernelPolyFit(dr);
    //return boundaryKernelChebFit(dr);
    return boundaryKernelAnalytic(dr);

}
std::pair<vec,scalar> closestPointTriangle(vec P, Triangle tri);
// returns the gradient for a flat infinite boundary at distance dr.
// Uses a simple numerical forward difference scheme
inline auto boundaryKernelDerivative(scalar d, scalar dx = (scalar)1e-7) {
    auto c = boundaryKernel(-d);
    auto p = boundaryKernel(-d + dx);
    return -(p - c) / dx;
}
std::tuple<bool, vec, scalar, scalar, vec> interactTriangle(vec p, Triangle tri);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentric(vec p, Triangle tri, scalar rhoi, scalar fi, scalar f0, scalar f1, scalar f2);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangles(vec p);
std::tuple<bool, vec, scalar, scalar, vec> interactLines(vec p, std::vector<vec> Poly, bool flipped = true);
std::tuple<bool, vec, scalar, scalar, vec> interactCircle(vec p);

// wrapper function around the boundary which will call the given function with:
// the distance to the boundary (scalar)
// the boundary integral (scalar)
// the gradient boundary integral (vec)
// the closest point on the boundary (vec)
// for any boundary b that is closer than 1.0 to the given position 
template <typename Func> auto boundaryFunc(const Particle& ptcl, Func&& c) {
    auto p = ptcl.pos;
    if (simulationMethod == boundaryMethod::analytical) {
        for (auto ti : ptcl.neighborTriangles) {
            const auto& tri = triangles[ti];
            auto [hit, pb, d, k, gk] = interactTriangle(p, tri);
            if (hit)
                c(pb, d, k, gk, true);
        }
    }
    if (simulationMethod == boundaryMethod::semi) {
        auto [hit, pb, d, k, gk] = interactTriangles(p);
        if (hit)
            c(pb, d, k, gk, true);
    }
    if (simulationMethod == boundaryMethod::sdf) {
        auto [hit, pb, d, k, gk] = interactLines(p, polygon);
        if (hit)
            c(pb, d, k, gk, true);
        for (auto obs : obstacles) {
            auto [hit, pb, d, k, gk] = interactLines(p, obs, false);
            if (hit)
                c(pb, d, k, gk, true);
        }
    }
}

template <typename Func> auto boundaryFunc(vec p, Func&& c) {

    if (simulationMethod != boundaryMethod::sdf) {
        //if (p.x() - domainMin.x() < (scalar)scale)
        //  c(vec(domainMin.x(), p.y()), p.x() - domainMin.x(), boundaryKernel(-(p.x() - domainMin.x() ) / scale),
        //    boundaryKernelDerivative((p.x() - domainMin.x()) / scale) * vec((scalar)1.0, (scalar)0.0) / scale, false);
        //if (domainMax.x() - p.x() < (scalar)scale)
        //  c(vec(domainMax.x(), p.y()), domainMax.x() - p.x(), boundaryKernel(-(domainMax.x() - p.x() ) / scale),
        //    boundaryKernelDerivative((domainMax.x() - p.x())/ scale) * vec((scalar)-1.0, (scalar)0.0) / scale, false);
        //if (p.y() - domainMin.y() < (scalar)scale)
        //  c(
        //      vec(p.x(), domainMin.y()), 
        //      p.y() - domainMin.y(), 
        //      boundaryKernel(-(p.y() - domainMin.y() ) / scale),
        //      boundaryKernelDerivative((p.y() - domainMin.y()) / scale) * vec((scalar)0.0, (scalar)1.0) / scale, false);
        //if (domainMax.y() - p.y() < (scalar)scale)
        //  c(vec(p.x(), domainMax.y()), domainMax.y() - p.y(), boundaryKernel(-(domainMax.y() - p.y() ) / scale),
        //    boundaryKernelDerivative((domainMax.y() - p.y()) / scale) * vec((scalar)0.0, (scalar)-1.0) / scale, false);
    }
    if (simulationMethod == boundaryMethod::analytical) {
        //auto [hit, pb, d, k, gk] = interactLines(p, polygon);

        //if (hit)
        //    c(pb, d, k, gk, true);

        for (auto& tri : triangles) {
            auto [hit, pb, d, k, gk] = interactTriangle(p, tri);
            if (hit)
                c(pb, d, k, gk, true);
        }
    }
    if (simulationMethod == boundaryMethod::semi) {
        auto [hit, pb, d, k, gk] = interactTriangles(p);
        if (hit)
            c(pb, d, k, gk, true);
    }
    if (simulationMethod == boundaryMethod::sdf) {

        //Polygon.push_back(P1);

        auto [hit, pb, d, k, gk] = interactLines(p, polygon);

        if (hit)
            c(pb, d, k, gk, true);
        for (auto obs : obstacles) {
            auto [hit, pb, d, k, gk] = interactLines(p, obs, false);

            if (hit)
                c(pb, d, k, gk, true);
        }
        {
            //auto [hit, pb, d, k, gk] = interactCircle(p);
            //if (hit)
            //    c(pb, d, k, gk, true);
        }
    }
}
