#include <cheb/cheb.h>
#include <cheb/algorithms.h>
#include <armadillo>
#include <numbers>
namespace cheb {
	inline vscalar fft(vscalar vals) {
		if (vals.size() <= 1)
			return vals;
		arma::vec v = arma::zeros(vals.size());
		for (int32_t i = 0; i < vals.size(); ++i)
			v[i] = vals[i];
		arma::cx_vec res = arma::fft(v);
		vscalar resv(res.size());
		for (int32_t i = 0; i < res.size(); ++i)
			resv[i] = res[i].real();
		return resv;
	}
	inline vscalar ifft(vscalar vals) {
		if (vals.size() <= 1)
			return vals;
		arma::cx_vec v(vals.size());
		for (int32_t i = 0; i < vals.size(); ++i)
			v[i] = vals[i];
		arma::cx_vec res = arma::ifft(v);
		vscalar resv(res.size());
		for (int32_t i = 0; i < res.size(); ++i)
			resv[i] = res[i].real();
		return resv;
	}
	vscalar rootsunit(vscalar ak, scalar htol) {
		auto n = standard_chop(ak);
		ak = subVec(ak, n);

		// Recurse for large coefficient sets
		if (n > 50) {
			auto chebpts = chebpts2((len)ak.size());
			auto lmap = Interval{ -1.,SPLITPOINT };
			auto rmap = Interval{ SPLITPOINT, 1. };
			auto lpts = lmap(chebpts);
			auto rpts = rmap(chebpts);
			auto lval = clenshaw(lpts, ak);
			auto rval = clenshaw(rpts, ak);
			auto lcfs = vals2coeffs2(lval);
			auto rcfs = vals2coeffs2(rval);
			auto lrts = rootsunit(lcfs, 2 * htol);
			auto rrts = rootsunit(rcfs, 2 * htol);
			auto ls = lmap(lrts);
			auto rs = rmap(rrts);
			return concatenate(lmap(lrts), rmap(rrts));
		}
		// A degree of 1 or less polynomial has no roots.
		// Note that this implicitly means that f(x) = 0 has no roots.
		if (n <= 1)
			return vscalar();
		vcomplex rts;
		// If the number of coefficients is 2 we can directly find the root
		// as the function described here is a linear function 
		if (n == 2)
			rts = vcomplex{ -ak[0] / ak[1] };
		// Otherweise determine the roots as the real eigenvalues of the
		// Chebyshev Companion Matrix
		else if (n <= 50) {
			arma::mat C = arma::zeros(n - 1, n - 1);
			C.diag(1) = arma::ones(n - 2) * .5;
			C.diag(-1) = arma::ones(n - 2) * .5;
			C(0, 1) = 1.;
			arma::mat D = arma::zeros(n - 1, n - 1);
			for (int32_t i = 0; i < n - 1; ++i)
				D(n - 2, i) = ak[i];
			arma::mat E = C - .5 * 1. / (*(std::end(ak) - 1)) * D;

			arma::cx_vec eigval;
			arma::cx_mat eigvec;
			arma::eig_gen(eigval, eigvec, E);

			rts = vcomplex(eigval.size());
			for (int32_t i = 0; i < eigval.size(); ++i)
				rts[i] = eigval(i);
		}
		// First we need to remove all roots with a complex part larger than some threshold
		auto mask = std::abs(evaluate(rts, [](auto c) {return c.imag(); })) < htol;
		auto rrts = evaluate(vcomplex(rts[mask]), [](auto c) {return c.real(); });
		// We then also remove any roots outside of the domain -1 1 (within some tolerance)
		// as these roots are not useful values and do not reflect actual roots.
		// Note that depending on the accuracy of the eigen solver, i.e. Eigen's eigen solver
		// this threshold needs to be fairly large to not miss any roots
		vscalar rrts2 = rrts[std::abs(rrts) <= 1. + 1. * htol];
		// Sort the roots from left to right for convienence
		rrts = sort(rrts2);
		// Ensure that the first and last root are in the correct interval to help the 
		// Newton polishing later
		if (rrts.size() >= 1) {
			rrts[0] = std::min(std::max(rrts[0], -1.), 1.);
			*(std::end(rrts) - 1) = std::min(*(std::end(rrts) - 1), 1.);
		}
		return rrts;
	}
	vscalar bary(vscalar xx, vscalar fk, vscalar xk, vscalar vk) {
		if (xx.size() == 0 || fk.size() == 0)
			return vscalar{};
		if (fk.size() == 1)
			return fk[0] * vscalar(1., xx.size());
		if (any(fk, [](auto s) {return std::isnan(s); }))
			return std::nan("") * vscalar(1., xx.size());
		vscalar out;
		// If the number of evaluation points is less than the number of nodes
		// we iterate over the evaluation points
		if (xx.size() < 4 * xk.size()) {
			out = vscalar(0.0, xx.size());
			for (int32_t i = 0; i < xx.size(); ++i) {
				auto tt = vk / (xx[i] - xk);
				out[i] = dot(tt, fk) / tt.sum();
			}
		}
		// else iterate over the nodes
		else {
			auto numer = vscalar(0.0, xx.size());
			auto denom = vscalar(0.0, xx.size());
			for (int32_t j = 0; j < xk.size(); ++j) {
				auto temp = vk[j] / (xx - xk[j]);
				numer = numer + temp * fk[j];
				denom = denom + temp;
			}
			out = numer / denom;
		}
		// In this process NaNs can arise, especially due to 
		// auto tt = vk / (xx[i] - xk);
		// These NaN's are always on function evaluation points that
		// coincide with the nodes of the polynomial. These NaNs can be resolved,
		// however, by setting the function values at these points to the function
		// values corresponding to the function values at the nodes (which can be
		// determined directly through the coefficients of the polynomial)
		for (int32_t k = 0; k < out.size(); ++k)
			if (std::isnan(out[k])) {
				auto idx = -1;
				for (int32_t i = 0; i < xk.size(); ++i)
					if (xx[k] == xk[i])
						idx = i;
				if (idx != -1)
					out[k] = fk[idx];
			}
		return out;
	}
	vscalar clenshaw(vscalar xx, vscalar ak) {
		if (xx.size() == 0 || ak.size() == 0)
			return vscalar{};
		if (ak.size() == 1)
			return ak[0] * vscalar(1., xx.size());
		if (any(ak, [](auto s) {return std::isnan(s); }))
			return std::nan("") * vscalar(1., xx.size());
		auto bk1 = 0 * xx;
		auto bk2 = 0 * xx;
		xx = 2. * xx;
		for (int32_t k = (len)ak.size() - 1; k > 1; k -= 2) {
			bk2 = ak[k] + xx * bk1 - bk2;
			bk1 = ak[k - 1] + xx * bk2 - bk1;
		}
		if ((ak.size() - 1) % 2 == 1) {
			auto temp = bk1;
			bk1 = ak[1] + xx * bk1 - bk2;
			bk2 = temp;
		}
		auto res = ak[0] + .5 * xx * bk1 - bk2;
		return res;
	}

	len round(scalar x)	{
		return (int)(x + 0.5);
	}

	len standard_chop(vscalar coeffs, scalar tol) {
		// This process only works if there are at least 17 coefficients
		auto n = (len)coeffs.size();
		auto cutoff = n;
		if (n < 17)
			return cutoff;
		// First convert the provided coefficients to a new vector envelope
		// that's monotically non-increasing (found using a scan of the flipped abs coefficients)
		auto b = flip(std::abs(coeffs));
		auto m = flip(scan(b, [](auto l, auto r) {return std::max(l, r); }));
		// If m[0] == 0. we cannot normalize, return 1 instead as fallback
		if (m[0] == 0.) return 1;
		// Normalize to begin at 0
		auto envelope = m / m[0];
		// Find the first plateau point, i.e. the point at which the coefficients are non-increasing
		auto plateauPoint = 0;
		int32_t j2 = 0;
		for (int32_t j = 1; j < n; ++j) {
			j2 = round(1.25 * (double)j + 5.);
			if (j2 > n - 1) // There is no plateau so all coefficients are required
				return cutoff;
			auto e1 = envelope[j];
			auto e2 = envelope[j2];
			auto r = 3 * (1. - std::log(e1) / std::log(tol));
			bool plateau = (e1 == 0.) || (e2 / e1 > r);
			if (plateau) { // found the first plateau point
				plateauPoint = j;
				break;
			}
		}
		//debugPrintVec(envelope);
		// The cutoff is the fixed to a point where the envelope, plus a linear function
		// included to bias the results towards the left end, is minimal.
		if (envelope[plateauPoint] == 0.)
			cutoff = plateauPoint;
		else {
			auto j3 = evaluate(envelope >= std::pow(tol, 7. / 6.), [](auto v) {return v ? 1 : 0; }).sum();
			if (j3 < j2) {
				j2 = j3 + 1;
				envelope[j2] = std::pow(tol, 7. / 6.);
			}
			auto cc = std::log10(subVec(envelope, j2));
			cc += linspace(0, (-1. / 3.) * std::log10(tol), j2);
			auto d = argmin(cc);
			cutoff = d;
		}
		return std::min(cutoff, n - 1);
	}
	vscalar adaptive(vfunc fun, uint32_t maxpow2) {
		vscalar coeffs; 
		for (int32_t k = 4; k < (len) maxpow2 + 1; ++k) {
			auto n = (len)std::pow(2, k) + 1;
			auto points = chebpts2(n);
			auto values = fun(points);
			coeffs = vals2coeffs2(values);

			auto chplen = standard_chop(coeffs);

			//auto midpoints = vscalar(0., n - 1);
			//for (int32_t i = 0; i < n - 1; ++i)
			//	midpoints[i] = (points[i] + points[i + 1]) * .5;
			//auto cvals = clenshaw(midpoints, coeffs);
			//auto fvals = fun(midpoints);
			//auto error = infnorm(fvals - cvals);
			//auto tol = 1e3 * eps;
			//if (error <= tol) {
			//	return slice(coeffs, 0, standard_chop(coeffs));
			//}

			if (chplen < coeffs.size())
				return slice(coeffs, 0, chplen);
			if (k == maxpow2) {
				std::cerr << "Constructor did not converge using " << n << " points" << std::endl;
				return coeffs;
			}
		}
		return coeffs;
	}
	vscalar newtonroots(cheb::IntervalFunction fun, vscalar rts, scalar tol, uint32_t maxiter) {
		auto [a, b] = fun.interval();
		if (rts.size() > 0) {
			auto dfun = fun.derivative() / (2. / (b - a));
			auto prv = std::numeric_limits<scalar>::infinity() * rts;
			auto count = 0;
			while (infnorm(rts - prv) > tol && count++ <= (len)maxiter) {
				prv = rts;
				rts = rts - fun.eval((rts)) / dfun.eval((rts));
			}
		}
		return rts;
	}
	vscalar coeffmult(vscalar fc, vscalar gc) {
		auto Fc = concatenate(2. * slice(fc, 0, 1), concatenate(slice(fc, 1), slice(fc, (len)fc.size() - 1, 0, -1)));
		auto Gc = concatenate(2. * slice(gc, 0, 1), concatenate(slice(gc, 1), slice(gc, (len)gc.size() - 1, 0, -1)));
		auto ak = ifft(fft(Fc) * fft(Gc));
		ak = concatenate(vscalar{ ak[0] }, slice(ak, 1) + slice(ak, (len)ak.size() - 1, 1, -1)) * 0.25;
		ak = subVec(ak, (len)fc.size());
		return ak;
	}
	vscalar barywts2(len n) {
		if (n == 0)
			return vscalar{};
		if (n == 1)
			return vscalar{ 1. };
		vscalar wts(1, n);
		*(std::end(wts) - 1) = .5;
		for (int32_t i = n - 2; i >= 0; i -= 2)
			wts[i] = -1;
		wts[0] = .5 * wts[0];
		return wts;
	}
	vscalar chebpts2(len n) {
		if (n == 1)
			return vscalar{ 0. };
		auto nn = flip(convert<vscalar>(arange(0, n)));
		auto pts = std::cos(nn * std::numbers::pi_v<scalar> / ((scalar)n - 1.));
		return pts;
	}
	vscalar vals2coeffs2(vscalar vals) {
		auto n = (len)vals.size();
		if (n <= 1) return vals;
		auto tmp = concatenate(flip(vals), slice(vals, 1, (len)vals.size() - 1));
		vscalar coeffs;
		coeffs = ifft(tmp);
		coeffs = slice(coeffs, 0, n);
		for (int32_t i = 1; i < n - 1; ++i) {
			coeffs[i] = 2. * coeffs[i];
		}
		return coeffs;
	}
	vscalar coeffs2vals2(vscalar coeffs) {
		auto n = (len)coeffs.size();
		if (n <= 1)
			return coeffs;
		for (int32_t i = 1; i < n - 1; ++i)
			coeffs[i] = .5 * coeffs[i];
		auto tmp = concatenate(coeffs, slice(coeffs, n - 2, 0, -1));
		vscalar values;
		values = fft(tmp);
		values = slice(values, n - 1, -1, -1);
		return values;
	}

}