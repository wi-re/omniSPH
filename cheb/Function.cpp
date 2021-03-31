#pragma once
#include <cheb/cheb.h>
#include <cheb/utilities.h>
#include <cheb/IntervalFunction.h>
#include <optional>
#include <variant>
namespace cheb {
	Function::Function(scalar c, Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(c, i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(std::string s, Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(vfunc f, std::variant<vlen, len> n, Domain d) {
		try {
			auto domain = d;
			auto nl = std::get<vlen>(n);
			auto is = domain.intervals();
			if (nl.size() == 1) {
				for (auto i : d.intervals())
					funs.emplace_back(f, nl[0], i);
				funs = check_funs(funs);
				breakdata = compute_breakdata(funs);
				return;
			}
			if (nl.size() != is.size())
				throw std::invalid_argument("Length of polynomial lengths and length of intervals don't match\n");
			for (int32_t i = 0; i < is.size(); ++i) {
				funs.push_back(IntervalFunction(f, nl[i], is[i]));
			}
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
		catch (std::bad_variant_access&) {
			auto nl = std::get<len>(n);
			for (auto i : d.intervals())
				funs.emplace_back(f, nl, i);
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
	}
	Function::Function(vfunc f, Domain d) {
		for (auto i : d.intervals())
			funs.emplace_back(f, i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(func fd, std::variant<vlen, len> n, Domain d) {
		auto f = funcToVfunc(fd);
		try {
			auto domain = d;
			auto nl = std::get<vlen>(n);
			auto is = domain.intervals();
			if (nl.size() == 1) {
				for (auto i : d.intervals())
					funs.emplace_back(f, nl[0], i);
				funs = check_funs(funs);
				breakdata = compute_breakdata(funs);
				return;
			}
			if (nl.size() != is.size())
				throw std::invalid_argument("Length of polynomial lengths and length of intervals don't match\n");
			for (int32_t i = 0; i < is.size(); ++i) {
				funs.push_back(IntervalFunction(f, nl[i], is[i]));
			}
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
		catch (std::bad_variant_access&) {
			auto nl = std::get<len>(n);
			for (auto i : d.intervals())
				funs.emplace_back(f, nl, i);
			funs = check_funs(funs);
			breakdata = compute_breakdata(funs);
		}
	}
	Function::Function(func fd, Domain d) {
		auto f = funcToVfunc(fd);
		for (auto i : d.intervals())
			funs.emplace_back(f, i);
		funs = check_funs(funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function(std::vector<IntervalFunction> _funs) {
		funs = check_funs(_funs);
		breakdata = compute_breakdata(funs);
	}
	Function::Function() {
		funs = check_funs({});
		breakdata = compute_breakdata(funs);
	}

	vscalar Function::operator()(vscalar x, method how) const {
		if (isempty()) return vscalar{};
		auto out = vscalar(std::nan(""), x.size());
		// Interior point values
		for (const auto& fun : funs) {
			auto idx = fun.interval().isinterior(x);
			auto fx = fun(x[idx], how);
			out[idx] = fx;
		}
		// Brekapoint values
		for (int32_t b = 0; b < breakdata.first.size(); ++b) 
			out[x == breakdata.first[b]] = breakdata.second[b];
		// Endpoint values
		auto lpts = x < breakdata.first[0];
		auto rpts = x > * (std::end(breakdata.first) - 1);
		out[lpts] = funs[0](x[lpts]);
		out[rpts] = (*std::rbegin(funs))(x[rpts]);
		return out;
	}
	scalar Function::operator()(scalar x, method how) const {
		if (isempty()) return 0.;
		return this->operator()(vscalar{ x }, how)[0];
	}
	bool Function::isempty() const {
		return funs.size() == 0;
	}
	vscalar Function::breakpoints()  const {
		return breakdata.first;
	}
	Domain Function::domain() const {
		if (isempty()) return Domain{ vscalar{} };
		return Domain(breakpoints());
	}
	Domain Function::support()const {
		if (isempty()) return Domain{ vscalar{} };
		return Domain{ {*std::begin(breakdata.first), *(std::end(breakdata.first) - 1)} };
	}
	scalar Function::hscale() const {
		if (isempty()) return 0.;
		return std::max(std::abs(*std::begin(breakdata.first)), std::abs(*(std::end(breakdata.first) - 1)));
	}
	bool Function::isconst() const {
		if (isempty()) return false;
		auto c = funs[0].coeffs()[0];
		for (const auto& fun : funs)
			if (!fun.isconst() || fun.coeffs()[0] != c) return false;
		return true;
	}
	scalar Function::vscale() const {
		if (isempty()) return 0.;
		scalar max = -DBL_MAX;
		for (const auto& fun : funs)
			max = std::max(max, fun.vscale());
		return max;
	}
	std::size_t Function::size() const {
		return funs.size();
	}
	Function Function::breakWith(Domain targetdomain) const {
		std::vector<IntervalFunction> newfuns;
		auto subintervals = targetdomain.intervals();
		auto intervalit = std::begin(subintervals);
		for (const auto& fun : funs)
			while (intervalit != std::end(subintervals) && fun.interval().contains(*intervalit))
				newfuns.push_back(fun.restricted(*intervalit++));
		return newfuns;
	}
	Function Function::simplified() const {
		if (isempty()) return Function();
		auto funs2 = funs;
		for (auto& fun : funs2)
			fun = fun.simplified();
		return Function(funs2);
	}
	Function Function::restricted(Domain subinterval, bool simplify) const {
		if (isempty()) return Function();
		auto newdom = domain().restricted(subinterval);
		Function fn = breakWith(newdom);;
		if (simplify) return fn.simplified();
		return fn;
	}
	vscalar Function::roots(bool polish) const {
		if (isempty()) return vscalar{};
		svec allrts;
		svec prvrts;
		auto htol = 1e2 * hscale() * eps;
		for (const auto& fun : funs) {
			auto rtsv = fun.roots(polish);
			auto rts = convert<svec>(rtsv);
			if (prvrts.size() > 0 && rts.size() > 0)
				if (std::abs(*std::rbegin(prvrts) - *std::begin(rts)) <= htol)
					rts = rts.size() == 1 ? svec{} : svec(std::begin(rts) + 1, std::end(rts));

			if (rts.size() > 0)
				allrts.insert(std::end(allrts), std::begin(rts), std::end(rts));
			prvrts = rts;
		}
		return convert<vscalar>(allrts);
	}
	Function Function::antiDerivative() const {
		std::vector<IntervalFunction> newfuns;
		std::optional<IntervalFunction> prevfun = std::nullopt;
		for (const auto& fun : funs) {
			auto integral = fun.antiDerivative();
			if (prevfun) {
				auto f = prevfun.value().endvalues();
				integral = integral + f[1];
			}
			newfuns.push_back(integral);
			prevfun = integral;
		}
		return Function(newfuns);
	}
	Function Function::derivative() const {
		std::vector<IntervalFunction> newfuns;
		for (const auto& fun : funs)
			newfuns.push_back(fun.derivative());
		return Function(newfuns);
	}
	scalar Function::definiteIntegral() const {
		scalar sum = 0.;
		for (const auto& fun : funs)
			sum += fun.defIntegral();
		return sum;
	}
	Function Function::operator-() {
		if (isempty()) return Function();
		std::vector<IntervalFunction> funs;
		for (auto& fun : this->funs)
			funs.push_back(-fun);
		return Function(funs);
	}
	scalar Function::dot(const Function& f) const {
		if (isempty()) return 0.;
		return (*this * f).definiteIntegral();
	}
	Function Function::absolute() {
		if (isempty()) return Function();
		auto newdom = domain().merged(roots());
		auto oldfuns = breakWith(newdom);
		std::vector<IntervalFunction> funs;
		for (auto& fun : oldfuns.funs)
			funs.push_back(fun.abs());
		return Function(funs);
	}
	Function Function::maximum(const Function& other) const {
		if (isempty() || other.isempty()) return Function();
		auto roots = (*this - other).roots();
		auto newdom = domain().united(other.domain()).merged(roots);
		auto interval = newdom.support();
		auto switch_ = Domain{ interval }.merged(roots);
		auto keys = .5 * (std::pow(-1., convert<vscalar>(arange(0, (len)switch_.size() - 1))) + 1);
		if (other(switch_[0]) > this->operator()(switch_[0]))
			keys = 1. - keys;
		std::vector<IntervalFunction> funs;
		auto ints = switch_.intervals();
		for (int32_t i = 0; i < ints.size(); ++i) {
			auto subdom = newdom.restricted(Domain{ ints[i].a, ints[i].b });
			Function subfun;
			if (keys[i])
				subfun = restricted(subdom);
			else
				subfun = other.restricted(subdom);
			funs.insert(std::end(funs), std::begin(subfun.funs), std::end(subfun.funs));
		}
		return Function(funs);
	}
	Function Function::maximum(const scalar& c) const {
		return maximum(Function(c, domain()));
	}
	Function Function::minimum(const Function& other) const {
		if (isempty() || other.isempty()) return Function();
		auto roots = (*this - other).roots();
		auto newdom = domain().united(other.domain()).merged(roots);
		auto interval = newdom.support();
		auto switch_ = Domain{ interval }.merged(roots);
		auto keys = .5 * (std::pow(-1., convert<vscalar>(arange(0, (len)switch_.size() - 1))) + 1);
		if (other(switch_[0]) < this->operator()(switch_[0]))
			keys = 1. - keys;
		std::vector<IntervalFunction> funs;
		auto ints = switch_.intervals();
		for (int32_t i = 0; i < ints.size(); ++i) {
			auto subdom = newdom.restricted(Domain{ ints[i].a, ints[i].b });
			Function subfun;
			if (keys[i])
				subfun = restricted(subdom);
			else
				subfun = other.restricted(subdom);
			funs.insert(std::end(funs), std::begin(subfun.funs), std::end(subfun.funs));
		}
		return Function(funs);
	}
	Function Function::minimum(const scalar& c) const {
		return minimum(Function(c, domain()));
	}
	Function pow(const Function& lhs, const Function& rhs) {
		if (lhs.isempty() || rhs.isempty())
			return Function();
		auto newdom = lhs.domain().united(rhs.domain());
		auto chbfn1 = lhs.breakWith(newdom);
		auto chbfn2 = rhs.breakWith(newdom);
		std::vector<IntervalFunction> newfuns;
		for (int32_t i = 0; i < chbfn1.funs.size(); ++i)
			newfuns.push_back(pow(chbfn1.funs[i], chbfn2.funs[i]).simplified());
		return Function(newfuns);
	}
	Function pow(const Function& lhs, const scalar& rhs) {
		if (lhs.isempty())
			return Function();
		auto chbfn1 = lhs;
		std::vector<IntervalFunction> newfuns;
		for (int32_t i = 0; i < chbfn1.funs.size(); ++i)
			newfuns.push_back(pow(chbfn1.funs[i], rhs));
		return Function(newfuns);
	}
	Function pow(const scalar& lhs, const Function& rhs) {
		if (rhs.isempty())
			return Function();
		auto chbfn1 = rhs;
		std::vector<IntervalFunction> newfuns;
		for (int32_t i = 0; i < chbfn1.funs.size(); ++i)
			newfuns.push_back(pow(lhs, chbfn1.funs[i]));
		return Function(newfuns);
	}
}