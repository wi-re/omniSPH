#include "2DMath.h"
#include <cfloat>
auto triangleArea(Triangle t) {
	auto [p0, p1, p2, i] = t; 
	return 0.5 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1]) + epsilon;
}
auto pointInTriangle(vec p, Triangle tri) {
	auto [p0, p1, p2, i] = tri;
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
}
vec closestPoint(vec P, vec A, vec B, bool clipped = false) {
	vec ap = P - A;
	vec ab = B - A;
	scalar ab2 = ab.dot(ab) + epsilon;
	scalar apab = ap.dot(ab);
	scalar t = apab / ab2;
	if (clipped)
		t = std::clamp(t, 0.0, 1.0);
	return A + ab * t;
}
auto closestPointEdge(vec c, vec p1, vec p2, vec center, bool check = true) {
	auto dC = (p2[1] - p1[1]) * c[0] - (p2[0] - p1[0]) * c[1] + p2[0] * p1[1] - p2[1] * p1[0];
	auto dT = (p2[1] - p1[1]) * center[0] - (p2[0] - p1[0]) * center[1] + p2[0] * p1[1] - p2[1] * p1[0];
	if (dC * dT < 0.0 || !check)
		return closestPoint(c, p1, p2);
	return c;
}
std::pair<vec, scalar> closestPointTriangle(vec P, Triangle tri) {
	if (pointInTriangle(P, tri)>=0) {
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
}

namespace trint {
using scalar = double;
using complex = std::complex<double>;
enum struct mode { scalar, gradient, diff, symmetric, basic };
enum struct variant {
  scalar,
  scalarBarycentric,
  gradient,
  gradientBarycentric,
  scalarGradient,
  scalarGradientBarycentric,
  all,
  allScalar
};
constexpr inline scalar epsilon((scalar)1e-11);

template <typename T = scalar>
struct vector2D {
  T x, y;
  auto norm() const { return std::sqrt(x * x + y * y); }
  auto operator*(const T& s) const { return vector2D<T>{x * s, y * s}; }
  auto operator/(const T& s) const { return vector2D<T>{x / s, y / s}; }
  auto operator-(const vector2D<T>& v) const {
    return vector2D<T>{x - v.x, y - v.y};
  }
  auto operator+(const vector2D<T>& v) const {
    return vector2D<T>{x + v.x, y + v.y};
  }
  auto dot(const vector2D<T>& r) const { return x * r.x + y * r.y; }
};
using vec2 = vector2D<>;
using vec2c = vector2D<complex>;

struct mat3 {
  scalar a00, a01, a02, a10, a11, a12, a20, a21, a22;
  template <typename T = scalar>
  auto mulVector(const vector2D<T>& v) const {
    return vector2D<T>{a00 * v.x + a01 * v.y, a10 * v.x + a11 * v.y};
  }
  template <typename T = scalar>
  auto mulPoint(const vector2D<T>& v) const {
    auto x = a00 * v.x + a01 * v.y + a02 * 1.0;
    auto y = a10 * v.x + a11 * v.y + a12 * 1.0;
    auto z = a20 * v.x + a21 * v.y + a22 * 1.0;
    return vector2D<T>{x / z, y / z};
  }
  inline mat3 transpose() const {
    return mat3{a00, a10, a20, a01, a11, a21, a02, a12, a22};
  }
  inline scalar det() const {
    return a00 * (a11 * a22 - a21 * a12) - a01 * (a10 * a22 - a12 * a20) +
           a02 * (a10 * a21 - a11 * a20);
  }
  inline mat3 invert() const {
    const auto invdet = 1 / det();
    const auto b00 = (a11 * a22 - a21 * a12) * invdet;
    const auto b01 = (a02 * a21 - a01 * a22) * invdet;
    const auto b02 = (a01 * a12 - a02 * a11) * invdet;
    const auto b10 = (a12 * a20 - a10 * a22) * invdet;
    const auto b11 = (a00 * a22 - a02 * a20) * invdet;
    const auto b12 = (a10 * a02 - a00 * a12) * invdet;
    const auto b20 = (a10 * a21 - a20 * a11) * invdet;
    const auto b21 = (a20 * a01 - a00 * a21) * invdet;
    const auto b22 = (a00 * a11 - a10 * a01) * invdet;
    return mat3{b00, b01, b02, b10, b11, b12, b20, b21, b22};
  }
};
struct Triangle {
  vec2 v0, v1, v2;
  scalar f0 = 1.0, f1 = 1.0, f2 = 1.0;
  scalar rho0 = 1.0, rho1 = 1.0, rho2 = 1.0;

  inline bool operator!=(const Triangle& rhs) const {
    return v0.x != rhs.v0.x || v0.y != rhs.v0.y || v1.x != rhs.v1.x ||
           v1.y != rhs.v1.y || v2.x != rhs.v2.x || v2.y != rhs.v2.y ||
           f0 != rhs.f0 || f1 != rhs.f1 || f2 != rhs.f2 || rho0 != rhs.rho0 ||
           rho1 != rhs.rho1 || rho2 != rhs.rho2;
  }

  inline Triangle transform(const mat3& T) const {
    return Triangle{T.mulPoint(v0),
                    T.mulPoint(v1),
                    T.mulPoint(v2),
                    f0,
                    f1,
                    f2,
                    rho0,
                    rho1,
                    rho2};
  }
  inline void print() const {
    std::cout << "Triangle with vertices:\n\t v0=[" << v0.x << "," << v0.y
              << "], v1=[" << v1.x << "," << v1.y << "], v2=[" << v2.x << ","
              << v2.y << "]\n";
    std::cout << "Quantities: f0=" << f0 << " f1=" << f1 << " f2=" << f2
              << "\n";
    std::cout << "Densities : r0=" << rho0 << " r1=" << rho1 << " r2=" << rho2
              << std::endl;
  }
  inline std::tuple<vec2, vec2, vec2> vertices() const {
    return std::make_tuple(v0, v1, v2);
  }
  inline std::tuple<scalar, scalar, scalar, scalar> baryWeights(
      const scalar x, const scalar y) const {
    const auto [x0, y0] = v0;
    const auto [x1, y1] = v1;
    const auto [x2, y2] = v2;
    const auto det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
    const auto l0 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / det;
    const auto l1 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / det;
    const auto l2 = 1. - l0 - l1;
    return std::make_tuple(l0, l1, l2, det);
  }
  inline scalar interp(const scalar x, const scalar y, const scalar fx0,
                       const scalar fx1, const scalar fx2) const {
    const auto [l0, l1, l2, d] = baryWeights(x, y);
    return l0 * fx0 + l1 * fx1 + l2 * fx2;
  }
  inline scalar interp(const scalar x, const scalar y) const {
    return interp(x, y, f0, f1, f2);
  }
  inline std::tuple<scalar, scalar, scalar, scalar> bigLambdas(
      const scalar x, const scalar y, const scalar fx0, const scalar fx1,
      const scalar fx2) const {
    const auto [x0, y0] = v0;
    const auto [x1, y1] = v1;
    const auto [x2, y2] = v2;
    const auto det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2);
    const auto l0 = (y1 - y2) * (fx0 - fx2) / det;
    const auto l1 = (x2 - x1) * (fx0 - fx2) / det;
    const auto l2 = (y2 - y0) * (fx1 - fx2) / det;
    const auto l3 = (x0 - x2) * (fx1 - fx2) / det;
    return std::make_tuple(l0, l1, l2, l3);
  }
  inline std::tuple<scalar, scalar, scalar, scalar> bigLambdas(
      const scalar x, const scalar y) const {
    return bigLambdas(x, y, f0, f1, f2);
  }
  inline std::tuple<scalar, scalar, scalar> factors(const scalar x,
                                                    const scalar y,
                                                    const scalar fx0,
                                                    const scalar fx1,
                                                    const scalar fx2) const {
    const auto [x2, y2] = v2;
    const auto [l0, l1, l2, l3] = bigLambdas(x, y, fx0, fx1, fx2);
    const auto t1 = l0 + l2;
    const auto t2 = l1 + l3;
    const auto t3 = -(l0 * x2 + l1 * y2 + l2 * x2 + l3 * y2) + fx2;
    if (std::abs(t1) > 1e-2 || std::abs(t2) > 1e-2) {
   //     std::cout << std::setprecision(4) << l0 << " : " << l1 << " : " << l2 << " : " << l3 << " -> " << t1 << " : " << t2 << " : " << t3 << " @ " << fx0 << " : " << fx1 << " : " << fx2 << " @ " << f0 << " : " << f1 << " : " << f2 << " @ " << rho0 << " : " << rho1 << " : " << rho2 << std::endl;
    }
    return std::make_tuple(t1, t2, t3);
  }
  inline std::tuple<scalar, scalar, scalar> factors(const scalar x,
                                                    const scalar y) const {
    return factors(x, y, f0, f1, f2);
  }

  inline std::tuple<scalar, scalar, scalar> taus(
      const scalar x, const scalar y, const mode m = mode::scalar,
      const scalar fi = 1.0, const scalar rhoi = 1.0) const {
    switch (m) {
      case mode::scalar:
        return factors(x, y);
      case mode::gradient:
        return factors(x, y);
      case mode::diff:
        return factors(x, y, rho0 * (f0 - fi) / rhoi, rho1 * (f1 - fi) / rhoi,
                       rho2 * (f2 - fi) / rhoi);
      case mode::symmetric:
        return factors(x, y,
                       rho0 * (f0 / (rho0 * rho0) + fi / (rhoi * rhoi)),
            rho1 * (f1 / (rho1 * rho1) + fi / (rhoi * rhoi)),
            rho2 * (f2 / (rho2 * rho2) + fi / (rhoi * rhoi)));
      case mode::basic:
        return std::make_tuple(0., 0., 1.);
    }
  }
};
struct integralState {
  Triangle t;

  scalar fi = 1.0;
  scalar rhoi = 1.0;

  mode m = mode::symmetric;

  std::tuple<scalar, scalar, scalar> sTaus;
  std::tuple<scalar, scalar, scalar> gTaus;

  integralState(const Triangle& _t, const scalar f, const scalar rho, mode _m)
      : t(_t), fi(f), rhoi(rho), m(_m) {
    sTaus = t.taus(0, 0, m == mode::basic ? m : mode::scalar, fi, rhoi);
    gTaus = t.taus(0, 0, m, fi, rhoi);
  }

  integralState transform(const mat3& T) const {
    auto TriT = t.transform(T);
    return integralState(TriT, fi, rhoi, m);
  }
};

constexpr inline auto sign(const scalar s) {
  if (s >= 0.) return 1.;
  return -1.;
}
template <int32_t n, typename T>
inline auto powers(const T v) {
  if constexpr (n == 0) return;
  if constexpr (n == 1) return v;
  const auto v2 = v * v;
  if constexpr (n == 2) return std::make_tuple(v, v2);
  const auto v3 = v2 * v;
  if constexpr (n == 3) return std::make_tuple(v, v2, v3);
  const auto v4 = v2 * v2;
  if constexpr (n == 4) return std::make_tuple(v, v2, v3, v4);
  const auto v5 = v3 * v2;
  if constexpr (n == 5) return std::make_tuple(v, v2, v3, v4, v5);
  const auto v6 = v3 * v3;
  if constexpr (n == 6) return std::make_tuple(v, v2, v3, v4, v5, v6);
  const auto v7 = v4 * v3;
  if constexpr (n == 7) return std::make_tuple(v, v2, v3, v4, v5, v6, v7);
  const auto v8 = v4 * v4;
  if constexpr (n == 8) return std::make_tuple(v, v2, v3, v4, v5, v6, v7, v8);
  const auto v9 = v5 * v4;
  if constexpr (n == 9)
    return std::make_tuple(v, v2, v3, v4, v5, v6, v7, v8, v9);
  const auto v10 = v5 * v5;
  if constexpr (n == 10)
    return std::make_tuple(v, v2, v3, v4, v5, v6, v7, v8, v9, v10);
  const auto v11 = v7 * v4;
  if constexpr (n == 11)
    return std::make_tuple(v, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11);
}
inline auto csqrt(const scalar val) {
  if (val >= 0.) return (complex)std::sqrt(val);
  return std::sqrt((complex)val);
}
inline auto csqrt(const complex val) {
  if (val.real() >= 0.) return (complex)std::sqrt(val.real());
  return std::sqrt(val);
}
template <typename T>
auto casin(const T val) {
  if (val >= -1. && val <= 1.) return std::asin(val);
  if (val < -1.) return -double_pi / 2.;
  if (val > 1.) return double_pi / 2.;
  return std::asin((complex)val).real();
}

inline auto lci(const vec2& origin, const vec2& direction, const vec2& center,
                const scalar r, const bool clipped = false,
                const scalar epsilon = DBL_EPSILON) {
  auto solveQuadratic = [epsilon, clipped](const scalar a, const scalar b,
                                           const scalar c) {
    const auto discr = b * b - 4. * a * c;
    auto x0 = 0., x1 = 0.;
    if (discr < 0)
      return std::make_tuple(false, x0, x1);
    else if (discr == 0.)
      x0 = x1 = -.5 * b / a;
    else {
      const auto q1 = -b + sqrt(discr);
      const auto q2 = -b - sqrt(discr);
      if (b >= 0) {
        x0 = q2 / (2. * a);
        x1 = (2. * c) / q2;
      } else {
        x1 = q1 / (2. * a);
        x0 = (2. * c) / q1;
      }
    }
    if (x0 > x1) {
      const auto t = x1;
      x1 = x0;
      x0 = t;
    }
    if (clipped)
      if (((-x0 > 1.0 + epsilon) && (-x1 > 1.0 + epsilon)) ||
          ((-x0 < 0.0 - epsilon) && (-x1 < 0.0 - epsilon)))
        return std::make_tuple(false, 0.0, 0.0);
    return std::make_tuple(true, x0, x1);
  };
  const auto [dx, dy] = direction;
  const auto [ox, oy] = origin;
  const auto [cx, cy] = center;
  const auto ocx = cx - ox;
  const auto ocy = cy - oy;
  const auto a = dx * dx + dy * dy;
  const auto b = 2.0 * (ocx * dx + ocy * dy);
  const auto c = (ocx * ocx + ocy * ocy) - r * r;
  const auto [hit, t0, t1] = solveQuadratic(a, b, c);
  if (!hit) return std::make_tuple(0, 0., 0.);
  return std::make_tuple(2, -t0, -t1);
};
scalar sgnEdge(const vec2& c, const vec2& p1, const vec2& p2) {
  return trint::sign((p2.y - p1.y) * c.x - (p2.x - p1.x) * c.y + p2.x * p1.y -
                     p2.y * p1.x);
}
scalar triangleArea(const std::tuple<vec2, vec2, vec2>& t) {
  const auto [p0, p1, p2] = t;
  return 0.5 * (-p1.y * p2.x + p0.y * (-p1.x + p2.x) + p0.x * (p1.y - p2.y) +
                p1.x * p2.y) +
         epsilon;
}
scalar triangleArea(const vec2& p0, const vec2& p1, const vec2& p2) {
  return 0.5 * (-p1.y * p2.x + p0.y * (-p1.x + p2.x) + p0.x * (p1.y - p2.y) +
                p1.x * p2.y) +
         epsilon;
}
scalar pointInTriangle(const vec2& p, const vec2& p0, const vec2& p1,
                       const vec2& p2) {
  const auto area = triangleArea(p0, p1, p2);
  const auto s =
      1.0 / (2.0 * area) *
      (p0.y * p2.x - p0.x * p2.y + (p2.y - p0.y) * p.x + (p0.x - p2.x) * p.y);
  const auto t =
      1.0 / (2.0 * area) *
      (p0.x * p1.y - p0.y * p1.x + (p0.y - p1.y) * p.x + (p1.x - p0.x) * p.y);
  const auto u = 1.0 - s - t;
  if (s > 0.0 && t > 0.0 && u > 0.0) {
    return 1.0;
  }
  if (s >= 0.0 && t >= 0.0 && u >= 0.0) {
    return 0.0;
  }
  return -1.0;
}
scalar pointInCircle(const vec2& c, const scalar r, const vec2& p) {
  const auto a = p.x - c.x;
  const auto b = p.y - c.y;
  const auto ab2 = a * a + b * b;
  const auto r2 = r * r;
  if (ab2 < r2) return 1.0;
  if (abs(ab2) / r2 < epsilon) return 0.0;
  return -1.0;
}
vec2 closestPoint(const vec2& P, const vec2& A, const vec2& B,
                  const bool clipped = false) {
  const vec2 ap = P - A;
  const vec2 ab = B - A;
  const scalar ab2 = ab.dot(ab) + epsilon;
  const scalar apab = ap.dot(ab);
  scalar t = apab / ab2;
  if (clipped) t = std::clamp(t, 0.0, 1.0);
  return A + ab * t;
}
namespace circle {
template <variant v = variant::scalarGradientBarycentric>
auto integrate(const integralState& state) {
  // These integral terms are specific to the kernel function being used and
  // need to be adjusted accordingly

  auto g1 = 0.;
  auto g2 = 0.;
  auto g3 = 1.;
  auto g1x = -1.;
  auto g2x = 0.;
  auto g3x = 0.;
  auto g1y = 0.;
  auto g2y = -1.;
  auto g3y = 0.;

  // End of Kernel function specific results

  {
    complex integral = 0.;
    vec2c gradient{0., 0.};
    // scalar, scalarBarycentric, gradient, gradientBarycentric, scalarGradient,
    // scalarGradientBarycentric, all
    if constexpr (v == variant::scalar || v == variant::scalarGradient)
      integral = g3;
    if constexpr (v == variant::scalarBarycentric ||
                  v == variant::scalarGradientBarycentric)
      integral = g1 * std::get<0>(state.sTaus) + g2 * std::get<1>(state.sTaus) +
                 g3 * std::get<2>(state.sTaus);
    if constexpr (v == variant::gradient || v == variant::scalarGradient)
      gradient = vec2c{g3x, g3y};
    if constexpr (v == variant::gradientBarycentric ||
                  v == variant::scalarGradientBarycentric)
      gradient = vec2c{
          g1x * std::get<0>(state.gTaus) + g2x * std::get<1>(state.gTaus) +
              g3x * std::get<2>(state.gTaus),
          g1y * std::get<0>(state.gTaus) + g2y * std::get<1>(state.gTaus) +
              g3y * std::get<2>(state.gTaus)};
    if constexpr (v == variant::scalar || v == variant::scalarBarycentric)
      return integral;
    if constexpr (v == variant::gradient || v == variant::gradientBarycentric)
      return gradient;
    if constexpr (v == variant::scalarGradient ||
                  v == variant::scalarGradientBarycentric)
      return std::make_pair(integral, gradient);
    if constexpr (v == variant::all)
      return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
    if constexpr (v == variant::allScalar)
      return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
  }
}
}  // namespace circle
namespace cone {
template <variant v = variant::scalarGradientBarycentric>
auto integrate(scalar cd, scalar gamma, const integralState& state) {
  // These integral terms are specific to the kernel function being used and
  // need to be adjusted accordingly
  const auto d = std::min(1.0, cd);
  auto [_d, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11] = powers<11>(d);
  auto c = cos(gamma);
  auto s = sin(gamma);
  auto s2 = sin(2. * gamma);

  auto g1 = (sin(gamma) * (525. * d11 - 3168. * d10 + 7700. * d9 - 9240. * d8 +
                           4950. * d7 - 924. * d5 + 165. * d3)) /
            (55. * double_pi);
  auto g2 = ((1. - c) * (525. * d11 - 3168. * d10 + 7700. * d9 - 9240. * d8 +
                         4950. * d7 - 924. * d5 + 165. * d3)) /
            (55. * double_pi);
  auto g3 = (gamma * (21. * d10 - 128. * d9 + 315. * d8 - 384. * d7 +
                      210. * d6 - 42. * d4 + 9. * d2)) /
            (2. * double_pi);
  auto g1x = ((s2 + 2. * gamma) * (84. * d10 - 448. * d9 + 945. * d8 -
                                   960. * d7 + 420. * d6 - 42. * d4)) /
             (4. * double_pi);
  auto g2x = ((1. - cos(gamma) * c) * (84. * d10 - 448. * d9 + 945. * d8 -
                                       960. * d7 + 420. * d6 - 42. * d4)) /
             (2. * double_pi);
  auto g3x = (s * (280. * d9 - 1512. * d8 + 3240. * d7 - 3360. * d6 +
                   1512. * d5 - 168. * d3)) /
             (3. * double_pi);
  auto g1y = g2x;
  auto g2y = -((s2 - 2. * gamma) * (84. * d10 - 448. * d9 + 945. * d8 -
                                    960. * d7 + 420. * d6 - 42. * d4)) /
             (4 * double_pi);
  auto g3y = ((1. - c) * (280. * d9 - 1512. * d8 + 3240. * d7 - 3360. * d6 +
                          1512. * d5 - 168. * d3)) /
             (3. * double_pi);

  // End of Kernel function specific results

  {
    complex integral = 0.;
    vec2c gradient{0., 0.};
    // scalar, scalarBarycentric, gradient, gradientBarycentric, scalarGradient,
    // scalarGradientBarycentric, all
    if constexpr (v == variant::scalar || v == variant::scalarGradient)
      integral = g3;
    if constexpr (v == variant::scalarBarycentric ||
                  v == variant::scalarGradientBarycentric)
      integral = g1 * std::get<0>(state.sTaus) + g2 * std::get<1>(state.sTaus) +
                 g3 * std::get<2>(state.sTaus);
    if constexpr (v == variant::gradient || v == variant::scalarGradient)
      gradient = vec2c{g3x, g3y};
    if constexpr (v == variant::gradientBarycentric ||
                  v == variant::scalarGradientBarycentric)
      gradient = vec2c{
          g1x * std::get<0>(state.gTaus) + g2x * std::get<1>(state.gTaus) +
              g3x * std::get<2>(state.gTaus),
          g1y * std::get<0>(state.gTaus) + g2y * std::get<1>(state.gTaus) +
              g3y * std::get<2>(state.gTaus)};
    if constexpr (v == variant::scalar || v == variant::scalarBarycentric)
      return integral;
    if constexpr (v == variant::gradient || v == variant::gradientBarycentric)
      return gradient;
    if constexpr (v == variant::scalarGradient ||
                  v == variant::scalarGradientBarycentric)
      return std::make_pair(integral, gradient);
    if constexpr (v == variant::all)
      return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
    if constexpr (v == variant::allScalar)
      return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
  }
}
}  // namespace cone
namespace planar {
template <variant v = variant::scalarGradientBarycentric>
auto integrate(scalar cd, const integralState& state) {
  // These integral terms are specific to the kernel function being used and
  // need to be adjusted accordingly
  const auto d = (scalar)std::clamp(cd, -1.0, 1.0);
  cd = d;
  const auto d2 = d * d;
  const auto d3 = d2 * d;
  const auto d4 = d2 * d2;
  const auto d5 = d3 * d2;
  const auto d6 = d3 * d3;
  const auto d7 = d4 * d3;
  const auto d8 = d4 * d4;
  const auto d9 = d5 * d4;
  const auto d10 = d5 * d5;
  const auto d11 = d7 * d4;
  auto ld = log((complex)d).real();
  auto l2d = log(2. * (complex)d).real();
  auto l_1 = log(complex(-1)).real();

  auto g1 = -((-10395. * d10 - 34650. * d8) * log(2. * sqrt(1. - d2) + 2.) +
              (10395. * d10 + 34650. * d8) * l2d +
              sqrt(1. - d2) * (2560. * d10 + 33125. * d8 + 12180. * d6 -
                               3556. * d4 + 832. * d2 - 96.)) /
            (330. * double_pi);
  auto g2 = 0.;
  auto g3 = abs(d) <= 1e-8 ? .5
                    : -((-525. * d9 - 1800. * d7) * log(sqrt(1. - d2) + 1.) +
                        (525. * d9 + 1800. * d7) * log(1. - sqrt(1. - d2)) -
                        30. * acos(d) +
                        sqrt(1. - d2) * (256. * d9 + 3398. * d7 + 1316. * d5 -
                                         420. * d3 + 130. * d)) /
                          (30. * double_pi);
  auto g1x = [=]() {
    if (abs(d) <= 1e-8) return (scalar)-.5;
    auto integral = scalar(0.0);
    if (d > 0.)
      return -((9450. * d9 + 25200. * d7) * log(abs(d)) +
               ((-26880. * d8 - 57600. * d6) * abs(d) + 5040. * d10 +
                56700. * d8 + 25200. * d6 - 2520. * d4) *
                   acos(d / abs(d)) +
               (-13125. * d9 - 34200. * d7) * log(sqrt(1. - d2) + 1.) +
               (3675. * d9 + 9000. * d7) *
                   log((complex)sqrt(1. - d2) - 1.).real() +
               60. * acos(d) +
               sqrt(1. - d2) * (4608. * d9 + 45984. * d7 + 11168. * d5 -
                                1880. * d3 + 60. * d) -
               3675. * l_1 * d9 - 9000. * l_1 * d7) /
             (60. * double_pi);
    auto d = (scalar)abs(cd);
    // const auto d2 = d * d;
    const auto d3 = d2 * d;
    // const auto d4 = d2 * d2;
    const auto d5 = d3 * d2;
    // const auto d6 = d3 * d3;
    const auto d7 = d4 * d3;
    // const auto d8 = d4 * d4;
    const auto d9 = d5 * d4;
    // const auto d10 = d5 * d5;
    const auto d11 = d7 * d4;
    // Note that the -1 is the result of circle::g1x()
    return -1. + ((9450. * d9 + 25200. * d7) * log(abs(d)) +
                  ((-26880. * d8 - 57600. * d6) * abs(d) + 5040. * d10 +
                   56700. * d8 + 25200. * d6 - 2520. * d4) *
                      acos(d / abs(d)) +
                  (-13125. * d9 - 34200. * d7) * log(sqrt(1. - d2) + 1.) +
                  (3675. * d9 + 9000. * d7) *
                      log((complex)sqrt(1. - d2) - 1.).real() +
                  60. * acos(d) +
                  sqrt(1. - d2) * (4608. * d9 + 45984. * d7 + 11168. * d5 -
                                   1880. * d3 + 60. * d) -
                  3675. * l_1 * d9 - 9000. * l_1 * d7) /
                     (60. * double_pi);
  }();
  auto g2x = 0.;
  auto g3x =
      abs(d) <= 1e-8 ? -16. / (3. * double_pi)
              : -(((-945. * d8 - 2520. * d6) * log(2. * sqrt(1. - d2) + 2.) +
                   (945. * d8 + 2520. * d6) * l2d +
                   sqrt(1. - d2) *
                       (256. * d8 + 2639. * d6 + 690. * d4 - 136. * d2 + 16.)) /
                  (3. * double_pi));
  auto g1y = g2x;
  auto g2y = abs(d) <= 1e-8 ? -.5
                     : ((-5775. * d9 - 16200. * d7) * log(sqrt(1. - d2) + 1.) +
                        (-3675. * d9 - 9000. * d7) *
                            log((complex)sqrt(1. - d2) - 1.).real() +
                        (9450. * d9 + 25200. * d7) * ld - 60. * acos(d) +
                        sqrt(1. - d2) * (512. * d9 + 6796. * d7 + 2632. * d5 -
                                         840. * d3 + 260. * d)) /
                           (60. * double_pi);
  auto g3y = 0.;

  //printf("\t%g -> [%g %g %g] [%g %g %g] [%g %g %g]\n", d, g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
  // End of Kernel function specific results

  {
    complex integral = 0.;
    vec2c gradient{0., 0.};
    // scalar, scalarBarycentric, gradient, gradientBarycentric, scalarGradient,
    // scalarGradientBarycentric, all
    if constexpr (v == variant::scalar || v == variant::scalarGradient)
      integral = g3;
    if constexpr (v == variant::scalarBarycentric ||
                  v == variant::scalarGradientBarycentric)
      integral = g1 * std::get<0>(state.sTaus) + g2 * std::get<1>(state.sTaus) +
                 g3 * std::get<2>(state.sTaus);
    if constexpr (v == variant::gradient || v == variant::scalarGradient)
      gradient = vec2c{g3x, g3y};
    if constexpr (v == variant::gradientBarycentric ||
                  v == variant::scalarGradientBarycentric)
      gradient = vec2c{
          g1x * std::get<0>(state.gTaus) + g2x * std::get<1>(state.gTaus) +
              g3x * std::get<2>(state.gTaus),
          g1y * std::get<0>(state.gTaus) + g2y * std::get<1>(state.gTaus) +
              g3y * std::get<2>(state.gTaus)};
    if constexpr (v == variant::scalar || v == variant::scalarBarycentric)
      return integral;
    if constexpr (v == variant::gradient || v == variant::gradientBarycentric)
      return gradient;
    if constexpr (v == variant::scalarGradient ||
                  v == variant::scalarGradientBarycentric)
      return std::make_pair(integral, gradient);
    if constexpr (v == variant::all)
      return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
    if constexpr (v == variant::allScalar)
      return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
  }
}
}  // namespace planar
namespace stub {
    template<variant v = variant::scalarGradientBarycentric>
    auto integrate(const scalar cd, const  scalar cl, const scalar beta, const scalar gamma, const integralState& state) {
        const auto d = (scalar)cd;
        const auto l = (scalar)cl;
        const auto s = sin(beta);
        const auto c = cos(beta);
        const auto u = (scalar)std::min(1., cl);
        const auto [d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11] = powers<11>(d);
        const auto [u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11] = powers<11>(u);
        const auto [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11] = powers<11>(l);
        const auto [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11] = powers<11>(s);
        const auto [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11] = powers<11>(c);
        auto [g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y] = cone::integrate<variant::allScalar>(d, gamma, state);

        auto csq_1i = csqrt(u2 - s2 * l2);
        auto csq_2i = csqrt(d2 - s2 * l2);

        auto l_1 = log(2. * csq_1i + 2. * u).real();
        auto l_2 = log(2. * csq_2i + 2. * d).real();
        auto l_3 = log((csq_1i + u) / u).real();
        auto l_4 = log(-(csq_1i - u) / u).real();
        auto l_5 = log((csq_2i + d) / d).real();
        auto l_6 = log(-(csq_2i - d) / d).real();

        auto csq_1 = csq_1i.real();
        auto csq_2 = csq_2i.real();

        auto cas_1 = casin((s * l) / u);
        auto cas_2 = casin((s * l) / d);

        g1 += (-((10395. * s11 * l10 + 34650. * s9 * l8) * l_1 + csq_1 * (6300. * s * u10 - 38016. * s * u9 + (92400. * s - 700. * s3 * l2) * u8 + (4752. * s3 * l2 - 110880. * s) * u7 + (-800. * s5 * l4 - 13200. * s3 * l2 + 59400. * s) * u6 + (5544. * s5 * l4 + 18480. * s3 * l2) * u5 + (-960. * s7 * l6 - 15840. * s5 * l4 - 11880. * s3 * l2 - 11088. * s) * u4 + (6930. * s7 * l6 + 23100. * s5 * l4) * u3 + (-1280. * s9 * l8 - 21120. * s7 * l6 - 15840. * s5 * l4 + 3696. * s3 * l2 + 1980. * s) * u2 + (10395. * s9 * l8 + 34650. * s7 * l6) * u - 2560. * s11 * l10 - 42240. * s9 * l8 - 31680. * s7 * l6 + 7392. * s5 * l4 - 1980. * s3 * l2) - 6930. * c * s * l * u10 + 42240. * c * s * l * u9 - 103950. * c * s * l * u8 + 126720. * c * s * l * u7 - 69300. * c * s * l * u6 + 13860. * c * s * l * u4 - 2970. * c * s * l * u2 + (-10395. * s11 * l10 - 34650. * s9 * l8) * l_2 + csq_2 * (2560. * s11 * l10 + (1280. * s9 * d2 - 10395. * s9 * d + 42240. * s9) * l8 + (960. * s7 * d4 - 6930. * s7 * d3 + 21120. * s7 * d2 - 34650. * s7 * d + 31680. * s7) * l6 + (800. * s5 * d6 - 5544. * s5 * d5 + 15840. * s5 * d4 - 23100. * s5 * d3 + 15840. * s5 * d2 - 7392. * s5) * l4 + (700. * s3 * d8 - 4752. * s3 * d7 + 13200. * s3 * d6 - 18480. * s3 * d5 + 11880. * s3 * d4 - 3696. * s3 * d2 + 1980. * s3) * l2 - 6300. * s * d10 + 38016. * s * d9 - 92400. * s * d8 + 110880. * s * d7 - 59400. * s * d6 + 11088. * s * d4 - 1980. * s * d2) + (6930. * c * s * d10 - 42240. * c * s * d9 + 103950. * c * s * d8 - 126720. * c * s * d7 + 69300. * c * s * d6 - 13860. * c * s * d4 + 2970. * c * s * d2) * l) / (660. * double_pi));
        g2 += (-((10395. * c * s10 * l10 + 34650. * c * s8 * l8) * l_1 - 6300. * u11 + csq_1 * (6300. * c * u10 - 38016. * c * u9 + (92400. * c - 700. * c * s2 * l2) * u8 + (4752. * c * s2 * l2 - 110880. * c) * u7 + (-800. * c * s4 * l4 - 13200. * c * s2 * l2 + 59400. * c) * u6 + (5544. * c * s4 * l4 + 18480. * c * s2 * l2) * u5 + (-960. * c * s6 * l6 - 15840. * c * s4 * l4 - 11880. * c * s2 * l2 - 11088. * c) * u4 + (6930. * c * s6 * l6 + 23100. * c * s4 * l4) * u3 + (-1280. * c * s8 * l8 - 21120. * c * s6 * l6 - 15840. * c * s4 * l4 + 3696. * c * s2 * l2 + 1980. * c) * u2 + (10395. * c * s8 * l8 + 34650. * c * s6 * l6) * u - 2560. * c * s10 * l10 - 42240. * c * s8 * l8 - 31680. * c * s6 * l6 + 7392. * c * s4 * l4 - 1980. * c * s2 * l2) + (6930. * s2 * l + 38016.) * u10 + (-42240. * s2 * l - 92400.) * u9 + (103950. * s2 * l + 110880.) * u8 + (-126720. * s2 * l - 59400.) * u7 + 69300. * s2 * l * u6 + 11088. * u5 - 13860. * s2 * l * u4 - 1980. * u3 + 2970. * s2 * l * u2 + (-10395. * c * s10 * l10 - 34650. * c * s8 * l8) * l_2 + csq_2 * (2560. * c * s10 * l10 + (1280. * c * s8 * d2 - 10395. * c * s8 * d + 42240. * c * s8) * l8 + (960. * c * s6 * d4 - 6930. * c * s6 * d3 + 21120. * c * s6 * d2 - 34650. * c * s6 * d + 31680. * c * s6) * l6 + (800. * c * s4 * d6 - 5544. * c * s4 * d5 + 15840. * c * s4 * d4 - 23100. * c * s4 * d3 + 15840. * c * s4 * d2 - 7392. * c * s4) * l4 + (700. * c * s2 * d8 - 4752. * c * s2 * d7 + 13200. * c * s2 * d6 - 18480. * c * s2 * d5 + 11880. * c * s2 * d4 - 3696. * c * s2 * d2 + 1980. * c * s2) * l2 - 6300. * c * d10 + 38016. * c * d9 - 92400. * c * d8 + 110880. * c * d7 - 59400. * c * d6 + 11088. * c * d4 - 1980. * c * d2) + (-6930. * s2 * d10 + 42240. * s2 * d9 - 103950. * s2 * d8 + 126720. * s2 * d7 - 69300. * s2 * d6 + 13860. * s2 * d4 - 2970. * s2 * d2) * l + 6300. * d11 - 38016. * d10 + 92400. * d9 - 110880. * d8 + 59400. * d7 - 11088. * d5 + 1980. * d3) / (660. * double_pi));
        g3 += (((-525. * s9 * l9 - 1800. * s7 * l7) * l_3 + (525. * s9 * l9 + 1800. * s7 * l7) * l_4 + (630. * cas_1 - 630. * beta) * u10 + (3840. * beta - 3840. * cas_1) * u9 + csq_1 * (70. * s * l * u8 - 480. * s * l * u7 + (80. * s3 * l3 + 1350. * s * l) * u6 + (-560. * s3 * l3 - 1920. * s * l) * u5 + (96. * s5 * l5 + 1620. * s3 * l3 + 1260. * s * l) * u4 + (-700. * s5 * l5 - 2400. * s3 * l3) * u3 + (128. * s7 * l7 + 2160. * s5 * l5 + 1680. * s3 * l3 - 420. * s * l) * u2 + (-1050. * s7 * l7 - 3600. * s5 * l5) * u + 256. * s9 * l9 + 4320. * s7 * l7 + 3360. * s5 * l5 - 840. * s3 * l3 + 270. * s * l) + (9450. * cas_1 - 9450. * beta) * u8 + (11520. * beta - 11520. * cas_1) * u7 + (6300. * cas_1 - 6300. * beta) * u6 + (1260. * beta - 1260. * cas_1) * u4 + (270. * cas_1 - 270. * beta) * u2 + (525. * s9 * l9 + 1800. * s7 * l7) * l_5 + (-525. * s9 * l9 - 1800. * s7 * l7) * l_6 + (-630. * d10 + 3840. * d9 - 9450. * d8 + 11520. * d7 - 6300. * d6 + 1260. * d4 - 270. * d2) * cas_2 + csq_2 * (-256. * s9 * l9 + (-128. * s7 * d2 + 1050. * s7 * d - 4320. * s7) * l7 + (-96. * s5 * d4 + 700. * s5 * d3 - 2160. * s5 * d2 + 3600. * s5 * d - 3360. * s5) * l5 + (-80. * s3 * d6 + 560. * s3 * d5 - 1620. * s3 * d4 + 2400. * s3 * d3 - 1680. * s3 * d2 + 840. * s3) * l3 + (-70. * s * d8 + 480. * s * d7 - 1350. * s * d6 + 1920. * s * d5 - 1260. * s * d4 + 420. * s * d2 - 270. * s) * l) + 630. * beta * d10 - 3840. * beta * d9 + 9450. * beta * d8 - 11520. * beta * d7 + 6300. * beta * d6 - 1260. * beta * d4 + 270. * beta * d2) / (60. * double_pi));
        g1x += (-(((4725. * s11 - 525. * s9) * l9 + (12600. * s9 - 1800. * s7) * l7) * l_3 + ((525. * s9 - 4725. * s11) * l9 + (1800. * s7 - 12600. * s9) * l7) * l_4 + (-2520. * cas_1 + 2520. * c * s + 2520. * beta) * u10 + (13440. * cas_1 - 13440. * c * s - 13440. * beta) * u9 + csq_1 * ((5600. * s3 - 3080. * s) * l * u8 + (16800. * s - 30240. * s3) * l * u7 + ((80. * s3 - 800. * s5) * l3 + (64800. * s3 - 36450. * s) * l) * u6 + ((5040. * s5 - 560. * s3) * l3 + (38400. * s - 67200. * s3) * l) * u5 + ((96. * s5 - 960. * s7) * l5 + (1620. * s3 - 12960. * s5) * l3 + (30240. * s3 - 17640. * s) * l) * u4 + ((6300. * s7 - 700. * s5) * l5 + (16800. * s5 - 2400. * s3) * l3) * u3 + ((128. * s7 - 1280. * s9) * l7 + (2160. * s5 - 17280. * s7) * l5 + (1680. * s3 - 10080. * s5) * l3 + (2100. * s - 3360. * s3) * l) * u2 + ((9450. * s9 - 1050. * s7) * l7 + (25200. * s7 - 3600. * s5) * l5) * u + (256. * s9 - 2560. * s11) * l9 + (4320. * s7 - 34560. * s9) * l7 + (3360. * s5 - 20160. * s7) * l5 + (3360. * s5 - 840. * s3) * l3) + (-28350. * cas_1 - 6300. * c * s3 * l2 + 28350. * c * s + 28350. * beta) * u8 + (28800. * cas_1 + 34560. * c * s3 * l2 - 28800. * c * s - 28800. * beta) * u7 + (-12600. * cas_1 - 75600. * c * s3 * l2 + 12600. * c * s + 12600. * beta) * u6 + 80640. * c * s3 * l2 * u5 + (1260. * cas_1 - 37800. * c * s3 * l2 - 1260. * c * s - 1260. * beta) * u4 + 5040. * c * s3 * l2 * u2 + ((525. * s9 - 4725. * s11) * l9 + (1800. * s7 - 12600. * s9) * l7) * l_5 + ((4725. * s11 - 525. * s9) * l9 + (12600. * s9 - 1800. * s7) * l7) * l_6 + (2520. * d10 - 13440. * d9 + 28350. * d8 - 28800. * d7 + 12600. * d6 - 1260. * d4) * cas_2 + csq_2 * ((2560. * s11 - 256. * s9) * l9 + ((1280. * s9 - 128. * s7) * d2 + (1050. * s7 - 9450. * s9) * d + 34560. * s9 - 4320. * s7) * l7 + ((960. * s7 - 96. * s5) * d4 + (700. * s5 - 6300. * s7) * d3 + (17280. * s7 - 2160. * s5) * d2 + (3600. * s5 - 25200. * s7) * d + 20160. * s7 - 3360. * s5) * l5 + ((800. * s5 - 80. * s3) * d6 + (560. * s3 - 5040. * s5) * d5 + (12960. * s5 - 1620. * s3) * d4 + (2400. * s3 - 16800. * s5) * d3 + (10080. * s5 - 1680. * s3) * d2 - 3360. * s5 + 840. * s3) * l3 + ((3080. * s - 5600. * s3) * d8 + (30240. * s3 - 16800. * s) * d7 + (36450. * s - 64800. * s3) * d6 + (67200. * s3 - 38400. * s) * d5 + (17640. * s - 30240. * s3) * d4 + (3360. * s3 - 2100. * s) * d2) * l) + (6300. * c * s3 * d8 - 34560. * c * s3 * d7 + 75600. * c * s3 * d6 - 80640. * c * s3 * d5 + 37800. * c * s3 * d4 - 5040. * c * s3 * d2) * l2 + (-2520. * c * s - 2520. * beta) * d10 + (13440. * c * s + 13440. * beta) * d9 + (-28350. * c * s - 28350. * beta) * d8 + (28800. * c * s + 28800. * beta) * d7 + (-12600. * c * s - 12600. * beta) * d6 + (1260. * c * s + 1260. * beta) * d4) / (60. * double_pi));
        g2x += (-((945. * c * s10 * l9 + 2520. * c * s8 * l7) * l_1 - 252. * s2 * u10 + 1344. * s2 * u9 + csq_1 * (560. * c * s2 * l * u8 - 3024. * c * s2 * l * u7 + (6480. * c * s2 * l - 80. * c * s4 * l3) * u6 + (504. * c * s4 * l3 - 6720. * c * s2 * l) * u5 + (-96. * c * s6 * l5 - 1296. * c * s4 * l3 + 3024. * c * s2 * l) * u4 + (630. * c * s6 * l5 + 1680. * c * s4 * l3) * u3 + (-128. * c * s8 * l7 - 1728. * c * s6 * l5 - 1008. * c * s4 * l3 - 336. * c * s2 * l) * u2 + (945. * c * s8 * l7 + 2520. * c * s6 * l5) * u - 256. * c * s10 * l9 - 3456. * c * s8 * l7 - 2016. * c * s6 * l5 + 336. * c * s4 * l3) + ((630. * s4 - 315. * s2) * l2 - 2835. * s2) * u8 + ((1728. * s2 - 3456. * s4) * l2 + 2880. * s2) * u7 + ((7560. * s4 - 3780. * s2) * l2 - 1260. * s2) * u6 + (4032. * s2 - 8064. * s4) * l2 * u5 + ((3780. * s4 - 1890. * s2) * l2 + 126. * s2) * u4 + (252. * s2 - 504. * s4) * l2 * u2 + (-945. * c * s10 * l9 - 2520. * c * s8 * l7) * l_2 + csq_2 * (256. * c * s10 * l9 + (128. * c * s8 * d2 - 945. * c * s8 * d + 3456. * c * s8) * l7 + (96. * c * s6 * d4 - 630. * c * s6 * d3 + 1728. * c * s6 * d2 - 2520. * c * s6 * d + 2016. * c * s6) * l5 + (80. * c * s4 * d6 - 504. * c * s4 * d5 + 1296. * c * s4 * d4 - 1680. * c * s4 * d3 + 1008. * c * s4 * d2 - 336. * c * s4) * l3 + (-560. * c * s2 * d8 + 3024. * c * s2 * d7 - 6480. * c * s2 * d6 + 6720. * c * s2 * d5 - 3024. * c * s2 * d4 + 336. * c * s2 * d2) * l) + ((315. * s2 - 630. * s4) * d8 + (3456. * s4 - 1728. * s2) * d7 + (3780. * s2 - 7560. * s4) * d6 + (8064. * s4 - 4032. * s2) * d5 + (1890. * s2 - 3780. * s4) * d4 + (504. * s4 - 252. * s2) * d2) * l2 + 252. * s2 * d10 - 1344. * s2 * d9 + 2835. * s2 * d8 - 2880. * s2 * d7 + 1260. * s2 * d6 - 126. * s2 * d4) / (6. * double_pi));
        g3x += (-((945. * s9 * l8 + 2520. * s7 * l6) * l_1 + csq_1 * (560. * s * u8 - 3024. * s * u7 + (6480. * s - 80. * s3 * l2) * u6 + (504. * s3 * l2 - 6720. * s) * u5 + (-96. * s5 * l4 - 1296. * s3 * l2 + 3024. * s) * u4 + (630. * s5 * l4 + 1680. * s3 * l2) * u3 + (-128. * s7 * l6 - 1728. * s5 * l4 - 1008. * s3 * l2 - 336. * s) * u2 + (945. * s7 * l6 + 2520. * s5 * l4) * u - 256. * s9 * l8 - 3456. * s7 * l6 - 2016. * s5 * l4 + 336. * s3 * l2) - 630. * c * s * l * u8 + 3456. * c * s * l * u7 - 7560. * c * s * l * u6 + 8064. * c * s * l * u5 - 3780. * c * s * l * u4 + 504. * c * s * l * u2 + (-945. * s9 * l8 - 2520. * s7 * l6) * l_2 + csq_2 * (256. * s9 * l8 + (128. * s7 * d2 - 945. * s7 * d + 3456. * s7) * l6 + (96. * s5 * d4 - 630. * s5 * d3 + 1728. * s5 * d2 - 2520. * s5 * d + 2016. * s5) * l4 + (80. * s3 * d6 - 504. * s3 * d5 + 1296. * s3 * d4 - 1680. * s3 * d3 + 1008. * s3 * d2 - 336. * s3) * l2 - 560. * s * d8 + 3024. * s * d7 - 6480. * s * d6 + 6720. * s * d5 - 3024. * s * d4 + 336. * s * d2) + (630. * c * s * d8 - 3456. * c * s * d7 + 7560. * c * s * d6 - 8064. * c * s * d5 + 3780. * c * s * d4 - 504. * c * s * d2) * l) / (6. * double_pi));
        g1y = g2x;
        g2y += ((((4725. * s11 - 4200. * s9) * l9 + (12600. * s9 - 10800. * s7) * l7) * l_3 + ((4200. * s9 - 4725. * s11) * l9 + (10800. * s7 - 12600. * s9) * l7) * l_4 + (2520. * cas_1 + 2520. * c * s - 2520. * beta) * u10 + (-13440. * cas_1 - 13440. * c * s + 13440. * beta) * u9 + csq_1 * ((5600. * s3 - 2520. * s) * l * u8 + (13440. * s - 30240. * s3) * l * u7 + ((720. * s3 - 800. * s5) * l3 + (64800. * s3 - 28350. * s) * l) * u6 + ((5040. * s5 - 4480. * s3) * l3 + (28800. * s - 67200. * s3) * l) * u5 + ((864. * s5 - 960. * s7) * l5 + (11340. * s3 - 12960. * s5) * l3 + (30240. * s3 - 12600. * s) * l) * u4 + ((6300. * s7 - 5600. * s5) * l5 + (16800. * s5 - 14400. * s3) * l3) * u3 + ((1152. * s7 - 1280. * s9) * l7 + (15120. * s5 - 17280. * s7) * l5 + (8400. * s3 - 10080. * s5) * l3 + (1260. * s - 3360. * s3) * l) * u2 + ((9450. * s9 - 8400. * s7) * l7 + (25200. * s7 - 21600. * s5) * l5) * u + (2304. * s9 - 2560. * s11) * l9 + (30240. * s7 - 34560. * s9) * l7 + (16800. * s5 - 20160. * s7) * l5 + (3360. * s5 - 2520. * s3) * l3) + (28350. * cas_1 - 6300. * c * s3 * l2 + 28350. * c * s - 28350. * beta) * u8 + (-28800. * cas_1 + 34560. * c * s3 * l2 - 28800. * c * s + 28800. * beta) * u7 + (12600. * cas_1 - 75600. * c * s3 * l2 + 12600. * c * s - 12600. * beta) * u6 + 80640. * c * s3 * l2 * u5 + (-1260. * cas_1 - 37800. * c * s3 * l2 - 1260. * c * s + 1260. * beta) * u4 + 5040. * c * s3 * l2 * u2 + ((4200. * s9 - 4725. * s11) * l9 + (10800. * s7 - 12600. * s9) * l7) * l_5 + ((4725. * s11 - 4200. * s9) * l9 + (12600. * s9 - 10800. * s7) * l7) * l_6 + (-2520. * d10 + 13440. * d9 - 28350. * d8 + 28800. * d7 - 12600. * d6 + 1260. * d4) * cas_2 + csq_2 * ((2560. * s11 - 2304. * s9) * l9 + ((1280. * s9 - 1152. * s7) * d2 + (8400. * s7 - 9450. * s9) * d + 34560. * s9 - 30240. * s7) * l7 + ((960. * s7 - 864. * s5) * d4 + (5600. * s5 - 6300. * s7) * d3 + (17280. * s7 - 15120. * s5) * d2 + (21600. * s5 - 25200. * s7) * d + 20160. * s7 - 16800. * s5) * l5 + ((800. * s5 - 720. * s3) * d6 + (4480. * s3 - 5040. * s5) * d5 + (12960. * s5 - 11340. * s3) * d4 + (14400. * s3 - 16800. * s5) * d3 + (10080. * s5 - 8400. * s3) * d2 - 3360. * s5 + 2520. * s3) * l3 + ((2520. * s - 5600. * s3) * d8 + (30240. * s3 - 13440. * s) * d7 + (28350. * s - 64800. * s3) * d6 + (67200. * s3 - 28800. * s) * d5 + (12600. * s - 30240. * s3) * d4 + (3360. * s3 - 1260. * s) * d2) * l) + (6300. * c * s3 * d8 - 34560. * c * s3 * d7 + 75600. * c * s3 * d6 - 80640. * c * s3 * d5 + 37800. * c * s3 * d4 - 5040. * c * s3 * d2) * l2 + (2520. * beta - 2520. * c * s) * d10 + (13440. * c * s - 13440. * beta) * d9 + (28350. * beta - 28350. * c * s) * d8 + (28800. * c * s - 28800. * beta) * d7 + (12600. * beta - 12600. * c * s) * d6 + (1260. * c * s - 1260. * beta) * d4) / (60. * double_pi));
        g3y += (-((945. * c * s8 * l8 + 2520. * c * s6 * l6) * l_1 - 560. * u9 + csq_1 * (560. * c * u8 - 3024. * c * u7 + (6480. * c - 80. * c * s2 * l2) * u6 + (504. * c * s2 * l2 - 6720. * c) * u5 + (-96. * c * s4 * l4 - 1296. * c * s2 * l2 + 3024. * c) * u4 + (630. * c * s4 * l4 + 1680. * c * s2 * l2) * u3 + (-128. * c * s6 * l6 - 1728. * c * s4 * l4 - 1008. * c * s2 * l2 - 336. * c) * u2 + (945. * c * s6 * l6 + 2520. * c * s4 * l4) * u - 256. * c * s8 * l8 - 3456. * c * s6 * l6 - 2016. * c * s4 * l4 + 336. * c * s2 * l2) + (630. * s2 * l + 3024.) * u8 + (-3456. * s2 * l - 6480.) * u7 + (7560. * s2 * l + 6720.) * u6 + (-8064. * s2 * l - 3024.) * u5 + 3780. * s2 * l * u4 + 336. * u3 - 504. * s2 * l * u2 + (-945. * c * s8 * l8 - 2520. * c * s6 * l6) * l_2 + csq_2 * (256. * c * s8 * l8 + (128. * c * s6 * d2 - 945. * c * s6 * d + 3456. * c * s6) * l6 + (96. * c * s4 * d4 - 630. * c * s4 * d3 + 1728. * c * s4 * d2 - 2520. * c * s4 * d + 2016. * c * s4) * l4 + (80. * c * s2 * d6 - 504. * c * s2 * d5 + 1296. * c * s2 * d4 - 1680. * c * s2 * d3 + 1008. * c * s2 * d2 - 336. * c * s2) * l2 - 560. * c * d8 + 3024. * c * d7 - 6480. * c * d6 + 6720. * c * d5 - 3024. * c * d4 + 336. * c * d2) + (-630. * s2 * d8 + 3456. * s2 * d7 - 7560. * s2 * d6 + 8064. * s2 * d5 - 3780. * s2 * d4 + 504. * s2 * d2) * l + 560. * d9 - 3024. * d8 + 6480. * d7 - 6720. * d6 + 3024. * d5 - 336. * d3) / (6. * double_pi));

        // End of Kernel function specific results

        {
            complex integral = 0.;
            vec2c gradient{ 0.,0. };
#define ARGS
            //scalar, scalarBarycentric, gradient, gradientBarycentric, scalarGradient, scalarGradientBarycentric, all
            if constexpr (v == variant::scalar || v == variant::scalarGradient)
                integral = g3;
            if constexpr (v == variant::scalarBarycentric || v == variant::scalarGradientBarycentric)
                integral = g1 * std::get<0>(state.sTaus) + g2 * std::get<1>(state.sTaus) + g3 * std::get<2>(state.sTaus);
            if constexpr (v == variant::gradient || v == variant::scalarGradient)
                gradient = vec2c{ g3x, g3y };
            if constexpr (v == variant::gradientBarycentric || v == variant::scalarGradientBarycentric)
                gradient = vec2c{
                    g1x * std::get<0>(state.gTaus) + g2x * std::get<1>(state.gTaus) + g3x * std::get<2>(state.gTaus),
                    g1y * std::get<0>(state.gTaus) + g2y * std::get<1>(state.gTaus) + g3y * std::get<2>(state.gTaus)
            };
            if constexpr (v == variant::scalar || v == variant::scalarBarycentric)
                return integral;
            if constexpr (v == variant::gradient || v == variant::gradientBarycentric)
                return gradient;
            if constexpr (v == variant::scalarGradient || v == variant::scalarGradientBarycentric)
                return std::make_pair(integral, gradient);
            if constexpr (v == variant::all)
                return std::make_tuple(g1, g2, g3, g1x, g2x, g3x, g1y, g2y, g3y);
#undef ARGS
        }
    }
}  // namespace stub

template <variant v = variant::scalarGradientBarycentric>
auto integrateCircle(const integralState& state) {
    return circle::integrate<v>(state);
}
template <variant v = variant::scalarGradientBarycentric>
auto integrateCone(vec2 v1, vec2 v2, scalar d, const integralState& state) {
    auto [x1, y1] = v1;
    auto alpha_1 = std::atan2(y1, x1);
    auto sgn = sign(v2.x * sin(-alpha_1) + v2.y * cos(-alpha_1));
    mat3 T{
               cos(alpha_1),       sin(alpha_1), 0.,
        -sgn * sin(alpha_1), sgn * cos(alpha_1), 0.,
                         0.,                 0., 1.};
    vec2 v0{0., 0.};
    vec2 vp0 = T.mulPoint(v0);
    vec2 vp1 = T.mulPoint(v1);
    vec2 vp2 = T.mulPoint(v2);
    auto a = vp1 / vp1.norm();
    auto b = vp2 / vp2.norm();
    auto eps = 1e-8;
    if (abs(a.x) < eps) a.x = 0.;
    if (abs(a.y) < eps) a.y = 0.;
    if (abs(b.x) < eps) b.x = 0.;
    if (abs(b.y) < eps) b.y = 0.;
    auto alpha = std::atan2(a.y, a.x);
    auto beta = std::atan2(b.y, b.x);
    auto gamma = beta - alpha;
    auto TriT = state.transform(T);
    auto [scalar, gradient] = cone::integrate<v>(d, gamma, TriT);
    auto Ti = T.invert();
    return std::make_pair(scalar, Ti.mulVector(gradient));
}
template <variant v = variant::scalarGradientBarycentric>
auto integratePlanar(vec2 v0, vec2 v1, vec2 v2, vec2 tri0, vec2 tri1,
                        const integralState& state) {
    vec2 p{0., 0.};
    auto cp_01 = closestPoint(p, tri0, tri1);
    auto d_01 = cp_01.norm();
    auto s_01c = sgnEdge(p, tri0, tri1);

    auto tri_c = (v0 + v1 + v2) / 3;
    auto s_01t = sgnEdge(tri_c, tri0, tri1);

    auto dist = complex(-s_01c * s_01t * d_01);

    auto [hits, l0, l1] = lci(tri0, tri1 - tri0, p, 1.);
    auto x0 = tri0 + (tri1 - tri0) * l0;
    auto x1 = tri0 + (tri1 - tri0) * l1;

    auto theta = acos(dist).real();
    auto theta1 = atan2(x0.y, x0.x);
    auto theta2 = atan2(x1.y, x1.x);
    auto flip = 1.;  // d_01 < 0. ? -1. : 1.;
    mat3 T{cos(theta - theta1),
            -sin(theta - theta1),
            0,
            flip * sin(theta - theta1),
            flip * cos(theta - theta1),
            0,
            0,
            0,
            1};
    auto xp1 = T.mulPoint(x1);
    auto xp0 = T.mulPoint(x0);
    auto s = abs(xp1.x - dist);
    if (abs(xp1.x - xp0.x) > abs(xp1.y - xp0.y)) {
    T = mat3{cos(theta - theta2),
                -sin(theta - theta2),
                0,
                flip * sin(theta - theta2),
                flip * cos(theta - theta2),
                0,
                0,
                0,
                1};
    xp1 = T.mulPoint(x1);
    xp0 = T.mulPoint(x1);
    }

    if (xp1.x * dist.real() < 0) {
    T = mat3{T.a00, T.a01, T.a02, -T.a10, -T.a11,
                T.a12, T.a20, T.a21, T.a22};
    xp1 = T.mulPoint(x1);
    xp0 = T.mulPoint(x0);
    }

    auto vp0 = T.mulPoint(v0);
    auto vp1 = T.mulPoint(v1);
    auto vp2 = T.mulPoint(v2);
    //printf("[%g %g] x [%g %g] -> [%g %g] x [%g %g] @ %g [%g %g %g]\n", x0.x, x0.y, x1.x, x1.y, xp0.x, xp0.y, xp1.x, xp1.y, dist.real(), theta, theta1, theta2);
    auto TriT = state.transform(T);
    auto Ti = T.invert();
    auto [scalar, gradient] = planar::integrate<v>(dist.real(), TriT);
    gradient = Ti.mulVector(gradient);

    return std::make_pair(scalar, gradient);
}
template <variant v = variant::scalarGradientBarycentric>
auto integrateTriangular(vec2 v1, vec2 v2, const integralState& state) {
    vec2 v0{0, 0};
    auto d1 = v1.norm();
    auto d2 = v2.norm();
    if (d1 < d2) {
    std::swap(v1, v2);
    }

    auto l = v1.norm();
    auto d = v2.norm();

    auto alpha = std::acos((v1 - v2).dot(v0 - v2) /
                            ((v1 - v2).norm() * (v0 - v2).norm()));
    auto beta = std::acos((v2 - v1).dot(v0 - v1) /
                        ((v0 - v1).norm() * (v2 - v1).norm()));
    auto gamma = std::acos((v1 - v0).dot(v2 - v0) /
                            ((v1 - v0).norm() * (v2 - v0).norm()));

    auto theta = std::atan2(v1.y, v1.x);
    auto s = sign(v2.x * sin(-theta) + v2.y * cos(-theta));

    mat3 T{cos(theta), sin(theta), 0, -s * sin(theta), s * cos(theta), 0,
            0,          0,          1};

    auto vp0 = T.mulPoint(v0);
    auto vp1 = T.mulPoint(v1);
    auto vp2 = T.mulPoint(v2);

    auto TriT = state.transform(T);

    auto [scalar, grad] =
        abs(alpha) >= epsilon && abs(beta) >= epsilon && abs(gamma) >= epsilon
            ? stub::integrate<v>(d, l, beta, gamma, TriT)
            : std::make_pair((complex)0., vec2c{0., 0.});

    auto gradient = T.invert().mulVector(grad);
    return std::make_pair(scalar, gradient);
}
template <variant v = variant::scalarGradientBarycentric>
auto integrateTri(vec2 v1, vec2 v2, const integralState& state) {
    vec2 v0{0, 0};
    auto d1 = v1.norm();
    auto d2 = v2.norm();
    if (d1 < d2) {
    std::swap(v1, v2);
    }
    auto alpha = std::acos((v1 - v2).dot(v0 - v2) /
                            ((v1 - v2).norm() * (v0 - v2).norm()));
    auto beta = std::acos((v2 - v1).dot(v0 - v1) /
                        ((v0 - v1).norm() * (v2 - v1).norm()));
    auto gamma = std::acos((v1 - v0).dot(v2 - v0) /
                            ((v1 - v0).norm() * (v2 - v0).norm()));

    if (gamma <= double_pi / 2. && alpha >= double_pi / 2. - 1e-7) {
    // Normal Case
    auto [i, g] = integrateTriangular<v>(v1, v2, state);
    return std::make_pair(i, g);

    } else {
    // Split Stub
    auto cp = closestPoint(v0, v1, v2);
    auto [i1, g1] = integrateTriangular<v>(v1, cp, state);
    auto [i2, g2] = integrateTriangular<v>(v2, cp, state);
    auto i = i1 + i2;
    auto g = vec2c{g1.x + g2.x, g1.y + g2.y};
    return std::make_pair(i, g);
    }
}
std::tuple<bool, scalar, scalar> intersectEdge(const vec2& q, const scalar r, const vec2& p, const vec2& d, const bool clip) {
    auto [hits, l1, l2] = lci(p, d, q, r, true, 1e-9);
    if (hits == 0) return std::make_tuple(false, 0., 0.);
    return std::make_tuple(true, l2, l1);

    const auto epsilon = 1e-9;
    const auto a = d.dot(d);
    const auto b = 2.0 * (d.dot(p) - d.dot(q));
    const auto c = p.dot(p) + q.dot(q) - 2.0 * p.dot(q) - r * r;
    if (b * b < 4.0 * a * c)
        return std::make_tuple(false, 0.0, 0.0);
    auto tn = (-b - std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a + epsilon);
    auto tp = (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a + epsilon);
    if (((tn > 1.0 + epsilon) && (tp > 1.0 + epsilon)) || ((tn < 0.0 - epsilon) && (tp < 0.0 - epsilon)))
        return std::make_tuple(false, 0.0, 0.0);
    if (clip) {
        tn = std::clamp(tn, 0.0, 1.0);
        tp = std::clamp(tp, 0.0, 1.0);
    }
    if (tn == tp)
        return std::make_tuple(false, 0.0, 0.0);
    return std::make_tuple(true, tn, tp);
}
template<variant v = variant::scalarGradientBarycentric>
auto integrateWedge(vec2 v0, vec2 v1, vec2 v2, const integralState& state) {
    vec2 p{ 0, 0 };
    if (triangleArea(v0, v1, v2) < 0.) {
        std::swap(v1, v2);
    }
    auto tri = std::make_tuple(v0, v1, v2);
    auto alpha = std::acos((v1 - v2).dot(v0 - v2) /
        ((v1 - v2).norm() * (v0 - v2).norm()));
    auto beta = std::acos((v2 - v1).dot(v0 - v1) /
        ((v0 - v1).norm() * (v2 - v1).norm()));
    auto gamma = std::acos((v1 - v0).dot(v2 - v0) /
        ((v1 - v0).norm() * (v2 - v0).norm()));
    auto l0 = v0.norm();
    auto l1 = v1.norm();
    auto l2 = v2.norm();
    //if (l0 < 1 || l1 < 1 || l2 < 1) {
    //    if (l2 < 1)
    //        tri = std::make_tuple(v2, v1, v0);
    //    else if (l1 < 1)
    //        tri = std::make_tuple(v1, v2, v0);
    //    else
    //        tri = std::make_tuple(v0, v1, v2);
    //}
    if (l0 < l1 && l0 < l2);
    else if (l1 < l0 && l1 < l2) {
        std::swap(v0, v1);
        tri = std::make_tuple(v0, v1, v2);
    }
    else if (l2 < l0 && l2 < l1) {
        std::swap(v0, v2);
        tri = std::make_tuple(v0, v1, v2);
    }
    else {
        if (alpha <= beta && alpha <= gamma)
            tri = std::make_tuple(v2, v1, v0);
        else if (beta <= alpha && beta <= gamma)
            tri = std::make_tuple(v1, v2, v0);
        else
            tri = std::make_tuple(v0, v1, v2);
    }
    if (triangleArea(tri) < 0.) {
        auto t = std::get<2>(tri);
        std::get<2>(tri) = std::get<1>(tri);
        std::get<1>(tri) = t;
    }
    auto [tri0, tri1, tri2] = tri;


    //auto [int1, tn1, tp1] = intersectEdge(p, 1., tri0, tri1 - tri0, false);
    //auto [int2, tn2, tp2] = intersectEdge(p, 1., tri0, tri2 - tri0, false);
    //auto [int3, tn3, tp3] = intersectEdge(p, 1., tri1, tri2 - tri1, false);

    auto [int1, tp1, tn1] = lci(tri0, tri1 - tri0, p, 1., true, 1e-9);
    auto [int2, tp2, tn2] = lci(tri0, tri2 - tri0, p, 1., true, 1e-9);
    auto [int3, tp3, tn3] = lci(tri1, tri2 - tri1, p, 1., true, 1e-9);

    if (int1 && !int2 && !int3)
        return integratePlanar<v>(tri0, tri1, tri2, tri1, tri0, state);
    if (!int1 && int2 && !int3)
        return integratePlanar<v>(tri0, tri1, tri2, tri0, tri2, state);
    if (!int1 && !int2 && int3)
        return integratePlanar<v>(tri0, tri1, tri2, tri2, tri1, state);

    auto t1 = tp1 <= 1. ? tp1 : tn1;
    auto t2 = tp2 <= 1. ? tp2 : tn2;
    auto t3 = tp3 <= 1. ? tp3 : tn3;

    auto xInter1 = tri0 + (tri1 - tri0) * t1;
    auto xInter2 = tri0 + (tri2 - tri0) * t2;
    auto xnInter1 = xInter1 / xInter1.norm();
    auto xnInter2 = xInter2 / xInter2.norm();

    auto aBase = triangleArea(p, xnInter1, xnInter2);
    auto sa = sign(aBase);
    auto pic = pointInTriangle(p, v0, v1, v2) >= 0.;

    complex integralCone = 0.;
    vec2c gradientCone{ 0., 0. };
    if (!int1 && !int2 && !int3 && pic) {
        auto [i, c] = integrateCircle<v>(state);
        integralCone = i;
        gradientCone = c;
    }
    else if (!int1 && !int2 && !int3 && !pic) {
    }
    else {
        auto [i, c] = integrateCone<v>(xnInter1, xnInter2, 1., state);
        integralCone = i;
        gradientCone = c;
    }
    auto a1 = triangleArea(p, tri0, xInter1);
    auto a2 = triangleArea(p, xInter2, tri0);
    auto [i1, g1] = int1 && ::abs(a1) >= 1e-7 ? integrateTri<v>(tri0, xInter1, state)
        : std::make_pair(complex(0.0), vec2c{ 0., 0. });
    auto [i2, g2] = int2 && ::abs(a2) >= 1e-7 ? integrateTri<v>(xInter2, tri0, state)
        : std::make_pair(complex(0.0), vec2c{ 0., 0. });
    auto integral = sa * integralCone + sign(a1) * i1 + sign(a2) * i2;
    auto gradient = gradientCone * sa + g1 * sign(a1) + g2 * sign(a2);

    auto [hits, lc1, lc2] = lci(tri1, tri2 - tri1, p, 1.);
    if ((hits > 0 && lc1 >= 0. && lc1 <= 1. && lc2 >= 0 && lc2 <= 1.) ||
        (hits == 0 && pic)) {
        auto [integralCircle, gradientCircle] = integrateCircle<v>(state);
        if (!int1 && !int2 && int3) {
            integralCircle = 0.;
            gradientCircle = vec2c{ 0., 0. };
        }
        if (!int1 && !int2 && !int3 && !pic) {
            return std::make_pair(complex(0.), vec2c{ 0., 0. });
        }
        if (!int1 && !int2 && !int3 && pic) {
            return std::make_pair(integralCircle, gradientCircle);
        }
        auto xc = tri1;
        xc.x = xc.x - (tri2 - tri1).y;
        xc.y = xc.y + (tri2 - tri1).x;
        auto s_01c = sgnEdge(xc, tri1, tri2);
        auto s_01p = sgnEdge(tri0, tri1, tri2);
        if (s_01c * s_01p < 0.) {
            xc.x = xc.x + 2. * (tri2 - tri1).y;
            xc.y = xc.y - 2. * (tri2 - tri1).x;
        }
        auto aPlanar = triangleArea(tri0, tri1, tri2);
        auto [integralPlanar, gradientPlanar] =
            integratePlanar<v>(xc, tri1, tri2, tri2, tri1, state);

        if (sa >= 0.) {
            integralPlanar = integralCircle - integralPlanar;
            gradientPlanar = gradientCircle - gradientPlanar;
        }
        integralCone = integralCone - sign(aPlanar) * integralPlanar;
        gradientCone = gradientCone - gradientPlanar * sign(aPlanar);
    }
    integral = sa * integralCone + sign(a1) * i1 + sign(a2) * i2;
    gradient = gradientCone * sa + g1 * sign(a1) + g2 * sign(a2);

    return std::make_pair(integral, gradient);
}
template<variant v = variant::scalarGradientBarycentric>
auto triIntegral(vec2 v0, vec2 v1, vec2 v2, const integralState& state) {
    auto triIntegral_0_1 = [&, state](vec2 v0, vec2 v1, vec2 v2) {
        if (triangleArea(v0, v1, v2) < 0.) {
            std::swap(v1, v2);
        }
        return integrateWedge<v>(v0, v1, v2, state);
    };
    auto triIntegral_2 = [&, state](vec2 v0, vec2 v1, vec2 v2) {
        if (triangleArea(v0, v1, v2) < 0.) {
            std::swap(v1, v2);
        }
        auto pic0 = pointInCircle(vec2{ 0, 0 }, 1., v0) >= 0;
        auto pic1 = pointInCircle(vec2{ 0, 0 }, 1., v1) >= 0;
        auto pic2 = pointInCircle(vec2{ 0, 0 }, 1., v2) >= 0;

        auto tri0 = std::make_tuple(v0, v1, v2);
        auto tri1 = std::make_tuple(v0, v1, v2);
        complex integral = 0.;
        vec2c gradient{ 0., 0. };
        if (pic0) {
            if (pic1) {
                auto d0 = (v1 - v0) / (v1 - v0).norm();
                auto vx1 = v0 + d0 * 3.;
                tri0 = std::make_tuple(v0, vx1, v2);
                tri1 = std::make_tuple(v1, vx1, v2);
            }
            if (pic2) {
                auto d0 = (v2 - v0) / (v2 - v0).norm();
                auto vx2 = v0 + d0 * 3.;
                tri0 = std::make_tuple(v0, v1, vx2);
                tri1 = std::make_tuple(v2, v1, vx2);
            }
        }
        else {
            auto d0 = (v2 - v1) / (v2 - v1).norm();
            auto vx1 = v1 + d0 * 3.;
            tri0 = std::make_tuple(v0, v1, vx1);
            tri1 = std::make_tuple(v0, v2, vx1);
        }

        auto [a0, a1, a2] = tri0;
        auto [b0, b1, b2] = tri1;
        auto [i0, g0] = triIntegral_0_1(a0, a1, a2);
        auto [i1, g1] = triIntegral_0_1(b0, b1, b2);
        integral = sign(triangleArea(a0, a1, a2)) * i0 - sign(triangleArea(b0, b1, b2)) * i1;
        gradient = g0 - g1;
        return std::make_pair(integral, gradient);
    };
    auto triIntegral_3 = [&, state](vec2 v0, vec2 v1, vec2 v2) {
        if (triangleArea(v0, v1, v2) < 0.)
            std::swap(v1, v2);
        auto l0 = (v0 - vec2{ 0,0 }).norm();
        auto l1 = (v1 - vec2{ 0,0 }).norm();
        auto l2 = (v2 - vec2{ 0,0 }).norm();

        if (l0 < l1 && l0 < l2);
        else if (l1 < l0 && l1 < l2) {
            std::swap(v0, v1);
            std::swap(v1, v2);
        }
        else if (l2 < l0 && l2 < l1) {
            std::swap(v0, v2);
            std::swap(v1, v2);
        }

        auto d1 = (v1 - v0) / (v1 - v0).norm();
        auto d2 = (v2 - v0) / (v2 - v0).norm();
        auto vx1 = v0 + d1 * 3.;
        auto vx2 = v0 + d2 * 3.;

        auto [i0, g0] = triIntegral_0_1(v0, vx1, vx2);
        auto [i1, g1] = triIntegral_2(v1, v2, vx1);
        auto [i2, g2] = triIntegral_0_1(v2, vx1, vx2);

        auto integral = i0 - i1 - i2;
        auto gradient = g0 - g1 - g2;

        return std::make_pair(integral, gradient);
    };

    auto c = 0;
    c += pointInCircle(vec2{ 0,0 }, 1., v0) >= 0 ? 1 : 0;
    c += pointInCircle(vec2{ 0,0 }, 1., v1) >= 0 ? 1 : 0;
    c += pointInCircle(vec2{ 0,0 }, 1., v2) >= 0 ? 1 : 0;

    if (c == 0 || c == 1) {
        return triIntegral_0_1(v0, v1, v2);
    }
    if (c == 2)
        return triIntegral_2(v0, v1, v2);
    return triIntegral_3(v0, v1, v2);

    return std::make_pair(complex(0), vec2c{ 0.,0. });
};
template<variant v = variant::scalarGradientBarycentric>
auto integrateTriangle(vec2 p, scalar r, vec2 v0, vec2 v1, vec2 v2, scalar fi, scalar ri, Triangle t, mode m) {
    auto l0 = (v0 - p).norm();
    auto l1 = (v1 - p).norm();
    auto l2 = (v2 - p).norm();
    if (l1 < l0 && l1 < l2) {
        auto t0 = v0;
        auto t1 = v1;
        auto t2 = v2;
        v0 = t1;
        v1 = t2;
        v2 = t0;
    }
    if (l2 < l0 && l2 < l1) {
        auto t0 = v0;
        auto t1 = v1;
        auto t2 = v2;
        v0 = t2;
        v1 = t0;
        v2 = t1;
    }


    if (triangleArea(v0, v1, v2) < 0.)
        std::swap(v1, v2);
    mat3 M{
        1. / r, 0.,-p.x / r,
        0.,1. / r, -p.y / r,
        0.,0.,1.
    };
    auto TriT = t.transform(M);
    auto vp0 = M.mulPoint(v0);
    auto vp1 = M.mulPoint(v1);
    auto vp2 = M.mulPoint(v2);
    auto pp = M.mulPoint(p);
    auto [integral, gradient] = triIntegral<v>(vp0, vp1, vp2, integralState(TriT, fi, ri, m));
    return std::make_pair(integral, gradient / r);
}

}

//
//std::tuple<bool, vec, scalar, scalar, vec> interactTriangleOld(vec p, Triangle tri) {
//    auto integral = std::clamp(integrateOverTriangle(p, 1.0, tri), 0.0, 1.0);
//    if (integral == 0.0) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
//    //std::cout << "Integral not equal 0!" << std::endl;
//    //std::cout << p << " -> " << integral << std::endl;
//
//    auto [pb, d] = closestPointTriangle(p, tri);
//    auto pentalty = 1.0;
//    scalar dx = 1e-3;
//    //std::cout << "Integral + dx" << std::endl;
//    auto integralxn = std::clamp(integrateOverTriangle(p - vec(dx, 0.0), 1.0, tri), 0.0, 1.0);
//    auto integralxp = std::clamp(integrateOverTriangle(p + vec(dx, 0.0), 1.0, tri), 0.0, 1.0);
//    //auto integralxn = integrateOverTriangle(p - vec(dx, 0.0), 1.0, tri);
//    //std::cout << "Integral + dy" << std::endl;
//    auto integralyn = std::clamp(integrateOverTriangle(p - vec(0.0, dx), 1.0, tri), 0.0, 1.0);
//    auto integralyp = std::clamp(integrateOverTriangle(p + vec(0.0, dx), 1.0, tri), 0.0, 1.0);
//    //auto integralyn = integrateOverTriangle(p - vec(0.0, dx), 1.0, tri);
//    vec gradient((integralxp - integralxn) / (2.0 * dx), (integralyp - integralyn) / (2.0 * dx));
//    //if(integral < 0.0 || integral > 1.0)
//    return std::make_tuple(true, pb, d, pentalty * integral, pentalty * gradient);
//}

std::tuple<bool, vec, scalar, scalar, vec> interactTriangle(vec p, scalar scale, Triangle tri) {
        auto [pb, d] = closestPointTriangle(p, tri);
        //if (d < -support) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));

        trint::vec2 v0{ tri.v0.x(), tri.v0.y() };
        trint::vec2 v1{ tri.v1.x(), tri.v1.y() };
        trint::vec2 v2{ tri.v2.x(), tri.v2.y() };
        auto [integral, gradient] = trint::integrateTriangle<trint::variant::scalarGradient>(trint::vec2{ p.x(),p.y() },  scale, v0, v1, v2, 1.0, 1.0, trint::Triangle{
            v0,v1,v2, 1.0,1.0,1.0, 1.0,1.0,1.0
            }, trint::mode::symmetric);
        vec grad = -vec(gradient.x.real(), gradient.y.real());
        if (integral == 0.0) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
       return std::make_tuple(true, pb, d, integral.real(), grad);
}
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentric(vec p, scalar scale, scalar rho0, Triangle tri, scalar rhoi, scalar fi, scalar f0, scalar f1, scalar f2) {
    auto [pb, d] = closestPointTriangle(p, tri);
    //if (d < -support) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));

    trint::vec2 v0{ tri.v0.x(), tri.v0.y() };
    trint::vec2 v1{ tri.v1.x(), tri.v1.y() };
    trint::vec2 v2{ tri.v2.x(), tri.v2.y() };
    auto [integral, gradient] = trint::integrateTriangle<trint::variant::scalarGradientBarycentric>(trint::vec2{ p.x(),p.y() }, scale, v0, v1, v2, fi, rhoi, trint::Triangle{
        v0,v1,v2, f0,f1,f2, rho0,rho0,rho0
        }, trint::mode::symmetric);
    vec grad = -vec(gradient.x.real(), gradient.y.real());
    if (integral == 0.0) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
    return std::make_tuple(true, pb, d, integral.real(), grad);
}
std::tuple<bool, vec, scalar, scalar, vec> interactTriangleBaryCentricSimple(vec p, scalar scale, scalar rho0, Triangle tri, scalar f0, scalar f1, scalar f2) {
    auto [pb, d] = closestPointTriangle(p, tri);
    //if (d < -support) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));

    trint::vec2 v0{ tri.v0.x(), tri.v0.y() };
    trint::vec2 v1{ tri.v1.x(), tri.v1.y() };
    trint::vec2 v2{ tri.v2.x(), tri.v2.y() };
    auto [integral, gradient] = trint::integrateTriangle<trint::variant::scalarGradientBarycentric>(trint::vec2{ p.x(),p.y() }, scale, v0, v1, v2, 0., 0., trint::Triangle{
        v0,v1,v2, f0,f1,f2, rho0,rho0,rho0
        }, trint::mode::scalar);
    vec grad = -vec(gradient.x.real(), gradient.y.real());
    if (integral == 0.0) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
    return std::make_tuple(true, pb, d, integral.real(), grad);
}

//
//std::tuple<bool, vec, scalar, scalar, vec> interactTriangles(vec p) {
//  vec pb_m;
//  scalar d_m = DBL_MAX;
//  for (auto &tri : triangles) {
//    auto [pb, d] = closestPointTriangle(p, tri);
//    // d = -(pb - p).norm();
//    if (pointInTriangle(p, tri))
//      d = -d;
//    //pb_m = pb;
//    //d_m = d;
//     if (d < d_m) {
//    	pb_m = pb;
//    	d_m = d;
//    }
//  }
//
//  if (d_m > 1.0)
//    return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
//  // std::cout << d_m << " -> " << p.x() << " . " << p.y() << " : " << pb_m.x() << " . " << pb_m.y() << std::endl;
//  vec x = p - pb_m;
//  auto dist = d_m;
//  if (x.norm() < 1e-5)
//    x = vec(0, 0);
//  else
//    x = x.normalized();
//
//  return std::make_tuple(true, pb_m, -dist, boundaryKernel(-dist),
//                         -boundaryKernelDerivative(-dist) * (dist < 0.f ? 1.f : -1.f) * x);
//}
//
//scalar sdpolygon(std::vector<vec> v, vec p) {
//  scalar d = (p - v[0]).dot(p - v[0]);
//  scalar s = 1.0;
//  for (int i = 0, j = v.size() - 1; i < v.size(); j = i, i++) {
//    vec e = v[j] - v[i];
//    vec w = p - v[i];
//    vec b = w - e * std::clamp(w.dot(e) / e.dot(e), 0.0, 1.0);
//    d = std::min(d, b.dot(b));
//    bool b1 = p.y() >= v[i].y();
//    bool b2 = p.y() < v[j].y();
//    bool b3 = e.x() * w.y() > e.y() * w.x();
//    if ((b1 && b2 && b3) || (!b1 && !b2 && !b3))
//      s *= -1.0;
//  }
//  return s * sqrt(d);
//}
//
//std::tuple<bool, vec, scalar, scalar, vec> interactLines(vec p, std::vector<vec> Poly, bool flipped) { 
//	scalar d = (flipped? -1.0 : 1.0) * sdpolygon(Poly, p);
//  scalar dx = (flipped ? -1.0 : 1.0) * sdpolygon(Poly, p + vec(0.001, 0.0));
//    scalar dy = (flipped ? -1.0 : 1.0) * sdpolygon(Poly, p + vec(0.0, 0.001));
//  vec grad = vec(dx - d, dy - d) / 0.001;
//    vec n = grad.normalized();
//
//	auto qd = d / support;
//    if (qd > 1.0)
//      return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
//    auto q = p - n * d;
//
//	auto gammaScale = 5.0;
//	auto gamma = (log(1.0 + exp(-gammaScale * qd)) - log(1.0 + exp(-gammaScale))) / log(2);
//	auto gammaPrime = -(gammaScale * exp(-gammaScale * qd)) / (log(2) * (exp(-gammaScale * qd) + 1.0));
//	gamma = 1.0;
//	gammaPrime = 0.0;
//
//	auto lambda = boundaryKernel(-qd);
//	auto lambdaPrime = boundaryKernelDerivative(-qd) * gamma + boundaryKernel(-qd) * gammaPrime;
//	lambda *= gamma;
//	//-(5 * e ^ (-(5 * x) / 2)) / (2 * ln(2) * (e ^ (-(5 * x) / 2) + 1))
//    //std::cout << d << " -> " << p.x() << " . " << p.y() << " : " << q.x() << " . " << q.y() << " : " <<grad.x() << " . " << grad.y() << " : " << n.x() << " . " << n.y() <<std::endl;
//    return std::make_tuple(true, q, -d, lambda,
//		lambdaPrime * n / support);
//}
//
//
//std::tuple<bool, vec, scalar, scalar, vec> interactCircle(vec p) {
//
//	auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
//	auto t = 12.5 / domainScale;
//
//	vec begin(domainEpsilon + 10 / domainScale, m - t);
//	vec end(domainWidth - 10 / domainScale, m + t); 
//	vec mid = (end + begin) * 0.5;
//	mid.x() -= 25.0 / domainScale;
//
//	auto thickness = 1.0 / domainScale;
//	auto c = 4. / domainScale;
//
//	vec center = mid;
//
//	vec dist = p - center;
//	scalar d = dist.norm() - 3.0 / domainScale;
//	auto grad = dist.normalized();
//	vec n = grad.normalized();
//
//	auto qd = d / support;
//	if (qd > 1.0)
//		return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
//	auto q = p - n * d;
//
//	auto gammaScale = 5.0;
//	auto gamma = (log(1.0 + exp(-gammaScale * qd)) - log(1.0 + exp(-gammaScale))) / log(2);
//	auto gammaPrime = -(gammaScale * exp(-gammaScale * qd)) / (log(2) * (exp(-gammaScale * qd) + 1.0));
//	gamma = 1.0;
//	gammaPrime = 0.0;
//
//	auto lambda = boundaryKernel(-qd);
//	auto lambdaPrime = boundaryKernelDerivative(-qd) * gamma + boundaryKernel(-qd) * gammaPrime;
//	lambda *= gamma;
//	//-(5 * e ^ (-(5 * x) / 2)) / (2 * ln(2) * (e ^ (-(5 * x) / 2) + 1))
//	//std::cout << d << " -> " << p.x() << " . " << p.y() << " : " << q.x() << " . " << q.y() << " : " <<grad.x() << " . " << grad.y() << " : " << n.x() << " . " << n.y() <<std::endl;
//	return std::make_tuple(true, q, -d, lambda,
//		lambdaPrime * n / support);
//
//}