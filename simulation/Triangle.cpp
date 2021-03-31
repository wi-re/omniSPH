#include "2DMath.h"
/*
# I:
# anti derivatives
# 0.0 -> 0.5: (4*a*(12*x^5-15*x^4+5*x^2))/(7*pi)
# 0.5 -> 1.0: -(4*a*x^2*(4*x^3-15*x^2+20*x-10))/(7*pi)

# definite integrals
# 0.0 -> l:  (4*a*(12*l^5-15*l^4+5*l^2))/(7*pi)

# 0 -> 0.5: (11*a)/(28*pi)
# 0.5 -> l: -(a*(64*l^5-240*l^4+320*l^3-160*l^2+13))/(28*pi)

# II:
# anti derivatives
# 0.0 -> 0.5: ((48*arcsin(sin(b)/x)-48*b)*x^5+(-18*sin(b)*(1-sin(b)^2/x^2)^(3/2)+30*sin(b)*sqrt(1-sin(b)^2/x^2)-60*arcsin(sin(b)/x)+60*b)*x^4-20*sin(b)*(1-sin(b)^2/x^2)^(3/2)*x^3+(20*arcsin(sin(b)/x)-20*b)*x^2+(20*sin(b)-60*sin(b)^3)*sqrt(1-sin(b)^2/x^2)*x+9*sin(b)^5*ln(sqrt(1-sin(b)^2/x^2)+1)-9*sin(b)^5*ln(sqrt(1-sin(b)^2/x^2)-1))/(7*pi)
# 0.5 -> 1.0: -((16*arcsin(sin(b)/x)-16*b)*x^5+(-6*sin(b)*(1-sin(b)^2/x^2)^(3/2)+10*sin(b)*sqrt(1-sin(b)^2/x^2)-60*arcsin(sin(b)/x)+60*b)*x^4+(-20*sin(b)*(1-sin(b)^2/x^2)^(3/2)+80*arcsin(sin(b)/x)-80*b)*x^3+(40*sin(b)*sqrt(1-sin(b)^2/x^2)-40*arcsin(sin(b)/x)+40*b)*x^2+(-60*sin(b)^3-40*sin(b))*sqrt(1-sin(b)^2/x^2)*x+(3*sin(b)^5+20*sin(b)^3)*ln(sqrt(1-sin(b)^2/x^2)+1)+(-3*sin(b)^5-20*sin(b)^3)*ln(sqrt(1-sin(b)^2/x^2)-1))/(7*pi)

# definite integrals
# l -> 0.5: (36*sin(b)^5*ln(abs(sqrt((l^2-sin(b)^2)/l^2)-1))-36*sin(b)^5*ln(sqrt((l^2-sin(b)^2)/l^2)+1)+(192*b-192*arcsin(sin(b)/l))*l^5+((l^2-sin(b)^2)/l^2)^(3/2)*(72*sin(b)*l^4+80*sin(b)*l^3)+sqrt((l^2-sin(b)^2)/l^2)*((240*sin(b)^3-80*sin(b))*l-120*sin(b)*l^4)+(240*arcsin(sin(b)/l)-240*b)*l^4+(80*b-80*arcsin(sin(b)/l))*l^2-36*sin(b)^5*ln(abs(sqrt(1-4*sin(b)^2)-1))+36*sin(b)^5*ln(sqrt(1-4*sin(b)^2)+1)+11*arcsin(2*sin(b))+sqrt(1-4*sin(b)^2)*(33*sin(b)-62*sin(b)^3)-11*b)/(28*pi)
# 0.5 -> 1: -((12*sin(b)^5+80*sin(b)^3)*ln(sqrt(1-sin(b)^2)+1)+(-12*sin(b)^5-80*sin(b)^3)*ln(1-sqrt(1-sin(b)^2))+(-12*sin(b)^5-80*sin(b)^3)*ln(sqrt(1-4*sin(b)^2)+1)+(12*sin(b)^5+80*sin(b)^3)*ln(1-sqrt(1-4*sin(b)^2))+13*arcsin(2*sin(b))-16*arcsin(sin(b))+sqrt(1-4*sin(b)^2)*(74*sin(b)^3+49*sin(b))+sqrt(1-sin(b)^2)*(-136*sin(b)^3-64*sin(b))+3*b)/(28*pi)

# l -> 1.0: ((3*sin(b)^5+20*sin(b)^3)*ln(sqrt((l^2-sin(b)^2)/l^2)+1)+(-3*sin(b)^5-20*sin(b)^3)*ln(sqrt((l^2-sin(b)^2)/l^2)-1)+(16*arcsin(sin(b)/l)-16*b)*l^5+sqrt((l^2-sin(b)^2)/l^2)*(10*sin(b)*l^4+40*sin(b)*l^2+(-60*sin(b)^3-40*sin(b))*l)+((l^2-sin(b)^2)/l^2)^(3/2)*(-6*sin(b)*l^4-20*sin(b)*l^3)+(60*b-60*arcsin(sin(b)/l))*l^4+(80*arcsin(sin(b)/l)-80*b)*l^3+(40*b-40*arcsin(sin(b)/l))*l^2+(-3*sin(b)^5-20*sin(b)^3)*ln(sqrt(1-sin(b)^2)+1)+(3*sin(b)^5+20*sin(b)^3)*ln(sqrt(1-sin(b)^2)-1)+4*arcsin(sin(b))+sqrt(1-sin(b)^2)*(34*sin(b)^3+16*sin(b))-4*b)/(7*pi)
*/

auto kernelAntiDerivative_0_05(scalar a, scalar x) {
	return (4.0 * a * (12.0 * power(x, 5) - 15.0 * power(x, 4) + 5.0 * power(x, 2))) / (7.0 * M_PI);
}
auto kernelAntiDerivative_05_1(scalar a, scalar x) {
	return  -(4.0 * a * power(x, 2) * (4.0 * power(x, 3) - 15.0 * power(x, 2) + 20.0 * x - 10.0)) / (7.0 * M_PI);
}
auto kernelIntegralS(scalar alpha, scalar l, scalar u) {
	if (l < 0.5)
		return (4.0 * alpha * (12.0 * power(l, 5) - 15.0 * power(l, 4) + 5 * power(l, 2))) / (7.0 * M_PI);
	else {
		auto a = (11.0 * alpha) / (28.0 * M_PI);
		auto b = -(alpha * (64.0 * power(l,5) - 240.0 * power(l,4) + 320.0 * power(l,3) - 160.0 * power(l,2) + 13.0)) / (28.0 * M_PI);
		return a + b;		
	}
}
auto angleAntiDerivative_0_05(scalar b, scalar xr, scalar d) {
	complex x(xr + epsilon, 0.0);

	auto a = (48.0 * std::asin(d * std::sin(b) / x) - 48.0 * b) * power(x, 5);

	auto b1 = -18.0 * d * std::sin(b) * std::pow(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2), 3.0 / 2.0);
	auto b2 = 30.0 * d * sin(b) * std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2));
	auto b3 = -60.0 * std::asin(d * std::sin(b) / x);
	auto b4 = 60.0 * b;

	auto c = 20.0 * d * std::sin(b) * std::pow(1.0 - power(d, 2) * power(sin(b), 2) / power(x, 2), 3.0 / 2.0);
	auto d1 = (20.0 * std::asin(d * std::sin(b) / x) - 20.0 * b);
	auto e = (20.0 * d * std::sin(b) - 60.0 * power(d, 3) * power(sin(b), 3)) * std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2));
	auto f = 9.0 * power(d, 5) * power(sin(b), 5) * std::log(std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2)) + 1.0);
	auto g = 9.0 * power(d, 5) * power(std::sin(b), 5) * std::log(std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2)) - 1.0);

	auto nom = (a + (b1 + b2 + b3 + b4) * power(x, 4) - c * power(x, 3) + d1 * power(x, 2) + e * x + f - g);
	auto denom = 7.0 * M_PI;
	return nom / denom;
}
auto angleAntiDerivative_05_1(scalar b, scalar xr, scalar d) {
	complex x(xr, 0.0);

	auto a = (16.0 * std::asin(d * sin(b) / x) - 16.0 * b) * power(x,5);

	auto b1 = -6.0 * d * std::sin(b) * std::pow(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2), 3.0 / 2.0);
	auto b2 = 10.0 * d * std::sin(b) * std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2));
	auto b3 = -60.0 * std::asin(d * std::sin(b) / x);
	auto b4 = 60.0 * b;

	auto c1 = -20.0 * d * std::sin(b) * std::pow(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2), 3.0 / 2.0) + 80.0 * std::asin(d * std::sin(b) / x);
	auto c2 = -80.0 * b;

	auto d1 = 40.0 * d * std::sin(b) * std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2));
	auto d2 = -40.0 * std::asin(d * std::sin(b) / x);
	auto d3 = 40.0 * b;

	auto e = (-60.0 * power(d,3) * power(std::sin(b),3) - 40.0 * d * std::sin(b)) * std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x,2));
	auto f = (3.0 * power(d,5) * power(std::sin(b),5) + 20.0 * power(d,3) * power(std::sin(b),3)) * std::log(std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2)) + 1.0);
	auto g = (-3.0 * power(d,5) * power(std::sin(b),5) - 20.0 * power(d,3)* power(std::sin(b),3)) * std::log(std::sqrt(1.0 - power(d, 2) * power(std::sin(b), 2) / power(x, 2)) - 1.0);

	auto nom = (a + (b1 + b2 + b3 + b4) * power(x, 4) + (c1 + c2) * power(x, 3) + (d1 + d2 + d3) * power(x, 2) + e * x + f + g);
	auto denom = -7.0 * M_PI;
	return nom / denom;
}
auto angleIntegralS(scalar b, scalar lower, scalar upper) {
	if (lower < 0.5) {
		if (upper < 0.5)
			return angleAntiDerivative_0_05(b, upper, upper) - angleAntiDerivative_0_05(b, lower, upper);
		auto a05u = angleAntiDerivative_0_05(b, 0.5, upper);
		auto a05l = angleAntiDerivative_0_05(b, lower, upper);
		auto a05 = a05u - a05l;
		auto a10u = angleAntiDerivative_05_1(b, upper, upper);
		auto a10l = angleAntiDerivative_05_1(b, 0.5, upper);
		auto a10 = a10u - a10l;
		return a05 + a10;
	}
	else {
      auto u = angleAntiDerivative_05_1(b, upper, upper);
      auto l = angleAntiDerivative_05_1(b, lower, upper);
      return u - l;
    }
}
auto integralSolution(scalar a, scalar b, scalar x, scalar l) {
	auto k = kernelIntegralS(a, x, l);
	auto ang = angleIntegralS(b, x, l);
	return k + ang;
}

auto angle(vec t0, vec t1, vec t2) {
	return std::acos((t1 - t0).dot(t2 - t0) / ((t1 - t0).norm() * (t2 - t0).norm()+epsilon));
}
auto beta(Triangle t) {
	auto [t0, t1, t2] = t;
	return std::acos((t1 - t0).dot(t2 - t0) / ((t1 - t0).norm() * (t2 - t0).norm() + epsilon));
}
auto sgnEdge(vec c, vec p1, vec p2) {
	return ::sgn((p2[1] - p1[1]) * c[0] - (p2[0] - p1[0]) * c[1] + p2[0] * p1[1] - p2[1] * p1[0]);
}
auto triangleArea(Triangle t) {
	auto [p0, p1, p2] = t; 
	return 0.5 * (-p1[1] * p2[0] + p0[1] * (-p1[0] + p2[0]) + p0[0] * (p1[1] - p2[1]) + p1[0] * p2[1]) + epsilon;
}
auto pointInTriangle(vec p, Triangle tri) {
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
}
auto pointInCircle(vec c, scalar r, vec p) {
	auto a = p.x() - c.x();
	auto b = p.y() - c.y();
    auto ab2 = a * a + b * b;
	auto r2 = r * r;
	if (ab2 < r2)
		return 1.0;
	if (abs(ab2) / r2 < epsilon)
		return 0.0;
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
auto closestPointTriangle(vec P, Triangle tri) {
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
auto intersectEdge(vec q, scalar r, vec p, vec d, bool clip = true) {
	auto epsilon = 1e-4;
	auto a = d.dot(d);
	auto b = 2.0 * (d.dot(p) - d.dot(q));
	auto c = p.dot(p) + q.dot(q) - 2.0 * p.dot(q) - r * r;
	if (b * b < 4.0 * a * c)
		return std::make_tuple(false, 0.0, 0.0);
    auto tn = (-b - std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a + epsilon);
	auto tp = (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a + epsilon);
	if(((tn > 1.0 + epsilon) && (tp > 1.0 + epsilon)) || ((tn <0.0 - epsilon) && (tp < 0.0 - epsilon)))
		return std::make_tuple(false, 0.0, 0.0);
	if (clip) {
		tn = std::clamp(tn, 0.0, 1.0);
		tp = std::clamp(tp, 0.0, 1.0);
	}
	if(tn == tp)
		return std::make_tuple(false, 0.0, 0.0);
	return std::make_tuple(true, tn, tp);
}

auto integrateTriangle(Triangle tri, vec p, bool Inverse = false, bool verbose = false) {
	if (triangleArea(tri) < 0.0)
		std::swap(tri.v1, tri.v2);
	auto [tri0, tri1, tri2] = tri;
	auto d = (tri0 - p).norm();
	if (d > 1.0) return 0.0;
	bool inside = false;
	if (pointInTriangle(p, tri) >= 0)
		inside = true;
	vec xInter1, xInter2;
	if (Inverse || inside) {
		auto [hit1, tn1, tp1] = intersectEdge(p, 1.0, tri0, tri0 - tri1, false);
		auto [hit2, tn2, tp2] = intersectEdge(p, 1.0, tri0, tri0 - tri2, false);
		xInter1 = tri0 + tp1 * (tri0 - tri1);
		xInter2 = tri0 + tp2 * (tri0 - tri2);
	}
	else {
		auto [hit1, tn1, tp1] = intersectEdge(p, 1.0, tri0, tri1 - tri0, false);
		auto [hit2, tn2, tp2] = intersectEdge(p, 1.0, tri0, tri2 - tri0, false);
		xInter1 = tri0 + tp1 * (tri1 - tri0);
		xInter2 = tri0 + tp2 * (tri2 - tri0);
	}
	auto alpha1 = angle(p, tri0, xInter1);
	auto alpha2 = angle(p, tri0, xInter2);
	auto beta1 = angle(xInter1, p, tri0);
	auto beta2 = angle(xInter2, p, tri0);
	
	auto a1 = triangleArea(Triangle{ p, tri0, xInter1 });
	auto a2 = triangleArea(Triangle{ p, xInter2, tri0 });

	complex integral1 = complex(0,0), integral2 = complex(0,0);
	if (std::abs(a2) >= 1e-6) {
		if ((alpha2 < M_PI / 2.0) && (M_PI - beta2 - alpha2 > 0.5 * M_PI - 1e-7)) {
			integral2 = integralSolution(alpha2, beta2, d, 1.0);
		}
		else {
			auto cp_12 = closestPoint(p, tri0, xInter2);
			auto d2 = (cp_12 - p).norm();
			auto alpha21 = angle(p, xInter2, cp_12);
			auto alpha22 = angle(p, cp_12, tri0);

			auto beta21 = angle(xInter2, p, cp_12);
			auto beta22 = angle(tri0, p, cp_12);

			auto integral21 = integralSolution(alpha21, beta21, d2, 1.0);
			auto integral22 = integralSolution(alpha22, beta22, d2, d);

			integral2 = integral21 + integral22;
		}
	}
	if (std::abs(a1) >= 1e-6) {
		if ((alpha1 < M_PI / 2.0) && (M_PI - beta1 - alpha1 > 0.5 * M_PI - 1e-7)) {
			integral1 = integralSolution(alpha1, beta1, d, 1.0);
		}
		else {
			auto cp_13 = closestPoint(p, tri0, xInter1);
			auto d3 = (cp_13 - p).norm();
			auto alpha11 = angle(p, xInter1, cp_13);
			auto alpha12 = angle(p, cp_13, tri0);

			auto beta11 = angle(xInter1, p, cp_13);
			auto beta12 = angle(tri0, p, cp_13);

			auto integral11 = integralSolution(alpha11, beta11, d3, 1.0);
			auto integral12 = integralSolution(alpha12, beta12, d3, d);

			integral1 = integral11 + integral12;
		}
	}
	auto alphaCone = angle(p, xInter1, xInter2);
	auto integralCone = alphaCone / (2.0 * M_PI);
	auto aBase = triangleArea(Triangle{ p,xInter1,xInter2 });
	auto integral = sgn(aBase) * integralCone + sgn(a1) * integral1 + sgn(a2) * integral2;

	auto cp_01 = closestPoint(p, tri0, tri1);
	auto cp_20 = closestPoint(p, tri0, tri2);
	
	auto d_01 = (cp_01 - p).norm();
	auto d_20 = (cp_20 - p).norm();

	auto k_01 = boundaryKernel(-d_01);
	auto k_20 = boundaryKernel(-d_20);

	if (inside)
		return 1.0 - k_01 - k_20 + integral.real();
	return std::abs(integral.real());
}
auto integrateOverTriangle(vec p, scalar r, Triangle t) {
	t.v0 /= scale;
	t.v1 /= scale;
	t.v2 /= scale;
	p /= scale;
	//r *= scale;
	if (triangleArea(t) < 0.0)
		std::swap(t.v1, t.v2);
	auto [tri0, tri1, tri2] = t;
	auto v0_visible = pointInCircle(p, r, tri0) >= 0;
	auto v1_visible = pointInCircle(p, r, tri1) >= 0;
	auto v2_visible = pointInCircle(p, r, tri2) >= 0;
	if (!v0_visible && !v1_visible && !v2_visible) 
	{
		auto [PC, dT] = closestPointTriangle(p, t);
		if (std::abs(dT) > r) {
			if (pointInTriangle(p, t) >= 0.0)
				return 1.0;
			return 0.0;
		}
		vec tri_c = (tri0 + tri1 + tri2) / 3.0;
		vec cp_01 = closestPoint(p, tri0, tri1);
		vec cp_12 = closestPoint(p, tri1, tri2);
		vec cp_20 = closestPoint(p, tri2, tri0);
		
		auto d_01 = (cp_01 - p).norm();
		auto d_12 = (cp_12 - p).norm();
		auto d_20 = (cp_20 - p).norm();

		auto A = boundaryKernel(-d_01);
		auto [hitA, pA, nA] = intersectEdge(p, r, tri0, tri1 - tri0);
		auto B = boundaryKernel(-d_12);
		auto [hitB, pB, nB] = intersectEdge(p, r, tri1, tri2 - tri1);
		auto C = boundaryKernel(-d_20);
		auto [hitC, pC, nC] = intersectEdge(p, r, tri2, tri0 - tri2);
		if (!hitA)
			A = 0.0;
		if (!hitB)
			B = 0.0;
		if (!hitC)
			C = 0.0;
		auto s_01c = (sgnEdge(p, tri0, tri1) >= 0.0 )? 1.0 : -1.0;
		auto s_12c = (sgnEdge(p, tri1, tri2) >= 0.0 )? 1.0 : -1.0;
		auto s_20c = (sgnEdge(p, tri2, tri0) >= 0.0 )? 1.0 : -1.0;
		auto s_01t = (sgnEdge(tri_c, tri0, tri1) >= 0.0) ? 1.0 : -1.0;
		auto s_12t = (sgnEdge(tri_c, tri1, tri2) >= 0.0) ? 1.0 : -1.0;
		auto s_20t = (sgnEdge(tri_c, tri2, tri0) >= 0.0) ? 1.0 : -1.0;
		//std::cout << p.x() << " : " << p.y() << "->" << s_01c << s_12c << s_20c << "." << s_01t << s_12t << s_20t << " -> " << A << " " << B << " " << C << ". " << (pointInTriangle(p, t) >= 0) << std::endl;
		if (pointInTriangle(p, t) >= 0)
			return 1.0 - A - B - C;
		return -s_01c * s_01t * A + -s_12c * s_12t * B + -s_20c * s_20t * C;
	}
	else 
	{
		//return 0.0;
		int32_t counter = 0;
		scalar integral = 0.0;
		if (v0_visible) {
			auto tr1 = tri0 + 20.0 * (tri1 - tri0);
			auto tr2 = tri0 + 20.0 * (tri2 - tri0);
			integral += integrateTriangle(Triangle{ tri0,tr1,tr2 }, p);
			counter++;
			if (!v2_visible && !v1_visible) {
				auto cp = closestPoint(p, tri1, tri2, true);
				auto d = cp.norm();
				if (d < 1.0)
					integral -= boundaryKernel(-d);
			}
		}
		if (v1_visible) {
			vec tr0 = tri1 + 20.0 * (tri0 - tri1);
			auto tr2 = tri1 + 20.0 * (tri2 - tri1);
			if (counter == 1)
				tr0 = tri1 - 20.0 * (tri0 - tri1);
			integral += (counter == 1 ? -1.0 : 1.0) * integrateTriangle(Triangle{ tri1,tr0,tr2 }, p);
			counter++;
			if (!v2_visible && !v0_visible) {
				auto cp = closestPoint(p, tri0, tri2, true);
				auto d = cp.norm();
				if (d < 1.0)
					integral -= boundaryKernel(-d);
			}
		}
		if (v2_visible) {
			vec tr0 = tri2 + 20.0 * (tri0 - tri2);
			vec tr1 = tri2 + 20.0 * (tri1 - tri2);
			if ((counter == 1 && v0_visible) || counter == 2)
				tr0 = tri2 - 20.0 * (tri0 - tri2);
			if ((counter == 1 && v1_visible) || counter == 2)
				tr1 = tri2 - 20.0 * (tri1 - tri2);
			integral += (counter == 1 ? -1.0 : 1.0) * integrateTriangle(Triangle{ tri2,tr0,tr1 }, p);
			counter++;
			if (!v0_visible && !v1_visible) {
				auto cp = closestPoint(p, tri1, tri0, true);
				auto d = cp.norm();
				if (d < 1.0)
					integral -= boundaryKernel(-d);
			}			
		}
		return integral;
	}
}

std::tuple<bool, vec, scalar, scalar, vec> interactTriangles2(vec p) {
	scalar dx = 1e-3;
	auto [hit, pb, d, integral, gk] = interactLines(p, polygon);
	auto [hit2, pb2, d2, integralxn, gk2] = interactLines(p - vec(dx, 0.0), polygon);
	auto [hit3, pb3, d3, integralxp, gk3] = interactLines(p + vec(dx, 0.0), polygon);
	auto [hit4, pb4, d4, integralyn, gk4] = interactLines(p - vec(0.0, dx), polygon);
	auto [hit5, pb5, d5, integralyp, gk5] = interactLines(p + vec(0.0, dx), polygon);
	for (auto tri : triangles) {

		auto integral2 = std::clamp(integrateOverTriangle(p, 1.0, tri), 0.0, 1.0);
		if (integral2 == 0.0) continue;
		integral += integral2;
		hit = true;
		//std::cout << "Integral not equal 0!" << std::endl;
		//std::cout << p << " -> " << integral << std::endl;

		auto [pb2, d2] = closestPointTriangle(p, tri);
		if (abs(d2) < abs(d)) {
			d = d2;
			pb = pb2;
		}
		auto pentalty = 1.0;
		//std::cout << "Integral + dx" << std::endl;
		integralxn += std::clamp(integrateOverTriangle(p - vec(dx, 0.0), 1.0, tri), 0.0, 1.0);
		integralxp += std::clamp(integrateOverTriangle(p + vec(dx, 0.0), 1.0, tri), 0.0, 1.0);
		//auto integralxn = integrateOverTriangle(p - vec(dx, 0.0), 1.0, tri);
		//std::cout << "Integral + dy" << std::endl;
		integralyn += std::clamp(integrateOverTriangle(p - vec(0.0, dx), 1.0, tri), 0.0, 1.0);
		integralyp += std::clamp(integrateOverTriangle(p + vec(0.0, dx), 1.0, tri), 0.0, 1.0);
		//auto integralyn = integrateOverTriangle(p - vec(0.0, dx), 1.0, tri);
	}
	//if(integral < 0.0 || integral > 1.0)
	//	std::cout << p.x() << " : " << p.y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << integral << ", " << integralxp << ", " << integralyp << " ==> " << gradient.x() << " : " << gradient.y() << std::endl;
	if (hit) {
		vec gradient((integralxp - integralxn) / (2.0 * dx), (integralyp - integralyn) / (2.0 * dx));
		return std::make_tuple(true, pb, d, integral, gradient);
	}
	return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
}

std::tuple<bool, vec, scalar, scalar, vec> interactTriangle(vec p, Triangle tri) {
	auto integral = std::clamp(integrateOverTriangle(p, 1.0, tri),0.0,1.0);
	if(integral == 0.0) return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
	//std::cout << "Integral not equal 0!" << std::endl;
	//std::cout << p << " -> " << integral << std::endl;
	
	auto [pb,d] = closestPointTriangle(p, tri);
	auto pentalty = 1.0;
	scalar dx = 1e-3;
	//std::cout << "Integral + dx" << std::endl;
	auto integralxn = std::clamp(integrateOverTriangle(p - vec(dx, 0.0), 1.0, tri), 0.0, 1.0);
	auto integralxp = std::clamp(integrateOverTriangle(p + vec(dx, 0.0), 1.0, tri),0.0,1.0);
	//auto integralxn = integrateOverTriangle(p - vec(dx, 0.0), 1.0, tri);
	//std::cout << "Integral + dy" << std::endl;
	auto integralyn = std::clamp(integrateOverTriangle(p - vec(0.0, dx), 1.0, tri), 0.0, 1.0);
	auto integralyp = std::clamp(integrateOverTriangle(p + vec(0.0, dx), 1.0, tri),0.0,1.0);
	//auto integralyn = integrateOverTriangle(p - vec(0.0, dx), 1.0, tri);
	vec gradient((integralxp - integralxn) / (2.0*dx), (integralyp - integralyn) / (2.0*dx));
	//if(integral < 0.0 || integral > 1.0)
	//	std::cout << p.x() << " : " << p.y() << " -> " << pb.x() << " : " << pb.y() << " @ " << d << " => " << integral << ", " << integralxp << ", " << integralyp << " ==> " << gradient.x() << " : " << gradient.y() << std::endl;

	return std::make_tuple(true, pb, d, pentalty * integral, pentalty * gradient );
}

std::tuple<bool, vec, scalar, scalar, vec> interactTriangles(vec p) {
  vec pb_m;
  scalar d_m = DBL_MAX;
  for (auto &tri : triangles) {
    auto [pb, d] = closestPointTriangle(p, tri);
    // d = -(pb - p).norm();
    if (pointInTriangle(p, tri))
      d = -d;
    //pb_m = pb;
    //d_m = d;
     if (d < d_m) {
    	pb_m = pb;
    	d_m = d;
    }
  }

  if (d_m > 1.0)
    return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
  // std::cout << d_m << " -> " << p.x() << " . " << p.y() << " : " << pb_m.x() << " . " << pb_m.y() << std::endl;
  vec x = p - pb_m;
  auto dist = d_m;
  if (x.norm() < 1e-5)
    x = vec(0, 0);
  else
    x = x.normalized();

  return std::make_tuple(true, pb_m, -dist, boundaryKernel(-dist),
                         -boundaryKernelDerivative(-dist) * (dist < 0.f ? 1.f : -1.f) * x);
}

scalar sdpolygon(std::vector<vec> v, vec p) {
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

std::tuple<bool, vec, scalar, scalar, vec> interactLines(vec p, std::vector<vec> Poly, bool flipped) { 
	scalar d = (flipped? -1.0 : 1.0) * sdpolygon(Poly, p);
  scalar dx = (flipped ? -1.0 : 1.0) * sdpolygon(Poly, p + vec(0.001, 0.0));
    scalar dy = (flipped ? -1.0 : 1.0) * sdpolygon(Poly, p + vec(0.0, 0.001));
  vec grad = vec(dx - d, dy - d) / 0.001;
    vec n = grad.normalized();

	auto qd = d / support;
    if (qd > 1.0)
      return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
    auto q = p - n * d;

	auto gammaScale = 5.0;
	auto gamma = (log(1.0 + exp(-gammaScale * qd)) - log(1.0 + exp(-gammaScale))) / log(2);
	auto gammaPrime = -(gammaScale * exp(-gammaScale * qd)) / (log(2) * (exp(-gammaScale * qd) + 1.0));
	gamma = 1.0;
	gammaPrime = 0.0;

	auto lambda = boundaryKernel(-qd);
	auto lambdaPrime = boundaryKernelDerivative(-qd) * gamma + boundaryKernel(-qd) * gammaPrime;
	lambda *= gamma;
	//-(5 * e ^ (-(5 * x) / 2)) / (2 * ln(2) * (e ^ (-(5 * x) / 2) + 1))
    //std::cout << d << " -> " << p.x() << " . " << p.y() << " : " << q.x() << " . " << q.y() << " : " <<grad.x() << " . " << grad.y() << " : " << n.x() << " . " << n.y() <<std::endl;
    return std::make_tuple(true, q, -d, lambda,
		lambdaPrime * n / support);
}


std::tuple<bool, vec, scalar, scalar, vec> interactCircle(vec p) {

	auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
	auto t = 12.5 / domainScale;

	vec begin(domainEpsilon + 10 / domainScale, m - t);
	vec end(domainWidth - 10 / domainScale, m + t); 
	vec mid = (end + begin) * 0.5;
	mid.x() -= 25.0 / domainScale;

	auto thickness = 1.0 / domainScale;
	auto c = 4. / domainScale;

	vec center = mid;

	vec dist = p - center;
	scalar d = dist.norm() - 3.0 / domainScale;
	auto grad = dist.normalized();
	vec n = grad.normalized();

	auto qd = d / support;
	if (qd > 1.0)
		return std::make_tuple(false, vec(0, 0), 0.0, 0.0, vec(0, 0));
	auto q = p - n * d;

	auto gammaScale = 5.0;
	auto gamma = (log(1.0 + exp(-gammaScale * qd)) - log(1.0 + exp(-gammaScale))) / log(2);
	auto gammaPrime = -(gammaScale * exp(-gammaScale * qd)) / (log(2) * (exp(-gammaScale * qd) + 1.0));
	gamma = 1.0;
	gammaPrime = 0.0;

	auto lambda = boundaryKernel(-qd);
	auto lambdaPrime = boundaryKernelDerivative(-qd) * gamma + boundaryKernel(-qd) * gammaPrime;
	lambda *= gamma;
	//-(5 * e ^ (-(5 * x) / 2)) / (2 * ln(2) * (e ^ (-(5 * x) / 2) + 1))
	//std::cout << d << " -> " << p.x() << " . " << p.y() << " : " << q.x() << " . " << q.y() << " : " <<grad.x() << " . " << grad.y() << " : " << n.x() << " . " << n.y() <<std::endl;
	return std::make_tuple(true, q, -d, lambda,
		lambdaPrime * n / support);

}