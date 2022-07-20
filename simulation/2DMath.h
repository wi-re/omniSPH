#pragma once
//#define _USE_MATH_DEFINES
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
inline scalar W(vec a, vec b, scalar ha, scalar hb) {
    vec r = b - a;
    auto support = (ha + hb)/2.;
    scalar q = r.norm() / support;
    //spline
    //if (q < scalar(0.5))
    //  return kernelConstant * (power(scalar(1.0) - q, 3) - scalar(4.0) * power(scalar(0.5) - q, 3));
    //else if (q <= scalar(1.0))
    //  return kernelConstant * power(scalar(1.0) - q, 3);
    if (q < 1.0)
        return 7.0 / double_pi / (support * support) * power(1.0 - q, 4) * (1.0 + 4.0 * q);
    // wend 4
    //if (q < 1.0)
    //    return 9.0 / pi / (support * support) * power(1.0 - q, 6) * (1.0 + 6.0 * q + 35.0 / 3.0 * q * q);
    return scalar(0.0);
}
// simple spiky kernel gradient function. Assumes H = 1.0
inline vec gradW(vec a, vec b, scalar ha, scalar hb) {
    vec r = a - b;
    auto support = (ha + hb)/2.;
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
    return -r / rn * 7.0 / double_pi / (support * support * support) * (20.0 * q * power(1 - q, 3));

    // wend4
    //return -r / rn * 9.0 / pi / (support * support * support) * 56.0 / 3.0 * q * (5.0 * q + 1.0) * power(1 - q, 5);

    // these lines represent a normal cubic spline gradient function, also assuming H = 1.0
    // if (q <= scalar(0.5))
    //  return r * kernelConstant * (-scalar(3.0) * power(scalar(1.0) - q, 2) + scalar(12.0) * power(scalar(0.5) - q, 2));
    // return -scalar(3.0) * r * kernelConstant * power(scalar(1.0) - q, 2);
}

std::pair<vec,scalar> closestPointTriangle(vec P, Triangle tri);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangle(vec p, scalar support, Triangle tri);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentric(vec p, scalar support, scalar rho0, Triangle tri, scalar rhoi, scalar fi, scalar f0, scalar f1, scalar f2);
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentricSimple(vec p, scalar support, scalar rho0, Triangle tri,scalar f0, scalar f1, scalar f2);