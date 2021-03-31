#include <cheb/cheb.h>
#include <cheb/vectorutil.h>
#include <random>
#include <sstream>

namespace cheb {
	bool all(const vbool& cfs) {
		bool init = true;
		for (const auto& val : cfs) {
			init = init && val;
		}
		return init;
	}
	bool all(const vscalar& cfs, std::function<bool(scalar)> predicate) {
		return all(evaluate(cfs, predicate));
	}
	bool all(const vcomplex& cfs, std::function<bool(complex)> predicate) {
		return all(evaluate(cfs, predicate));
	}
	bool any(const vbool& cfs) {
		bool init = false;
		for (const auto& val : cfs) {
			init = init || val;
		}
		return init;
	}
	bool any(const vscalar& cfs, std::function<bool(scalar)> func) {
		return any(evaluate(cfs, func));
	}
	bool any(const vcomplex& cfs, std::function<bool(complex)> func) {
		return all(evaluate(cfs, func));
	}
	scalar infnorm(const vscalar& vals) {
		return std::abs(vals).max();
	}
	vscalar randn(int32_t n) {
		static std::random_device rd;
		static std::default_random_engine gen(rd());
		static std::uniform_real_distribution<scalar> dis(0.0, 1.0);
		vscalar randData(n);
		for (int32_t i = 0; i < n; ++i)
			randData[i] = (dis(gen));
		return randData;
	}
	vscalar linspace(scalar min, scalar max, int32_t n) {
		vscalar data(n);
		for (int32_t i = 0; i < n; ++i)
			data[i] = (min + (max - min) / ((scalar)n - 1.0) * (scalar)i);
		return data;
	}
	vlen arange(len min, len max) {
		vlen valA(max - min);
		for (int32_t i = 0; i < max - min; ++i)
			valA[i] = min + i;
		return valA;
	}
	vfunc funcToVfunc(func fn) {
		return [fn](vscalar x) -> vscalar {return evaluate(x, fn); };
	}
}
