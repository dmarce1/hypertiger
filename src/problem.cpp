/*
 * problem.cpp
 *
 *  Created on: May 29, 2015
 *      Author: dmarce1
 */

#include "defs.hpp"
#include "problem.hpp"
#include "options.hpp"
#include "grid.hpp"
#include <cmath>
#include <utility>
#include <hpx/include/lcos.hpp>

constexpr integer euler_type::NF;
static constexpr integer D_i = 0;
static constexpr integer tau_i = 1;
static constexpr integer Sx_i = 2;
static constexpr integer Sy_i = 3;
static constexpr integer Sz_i = 4;

static constexpr integer rho_i = 0;
static constexpr integer eps_i = 1;
static constexpr integer ux_i = 2;
static constexpr integer uy_i = 3;
static constexpr integer uz_i = 4;

static const simd_vector zero(0.0);
static const simd_vector one(1.0);
static const simd_vector two(2.0);
static const simd_vector half(0.5);
static const simd_vector fgamma(7.0 / 5.0);

euler_type::vector_type euler_type::initial_value(
		const std::vector<simd_vector>& x, real) {
	std::vector<simd_vector> u(NF);
	for (integer i = 0; i != simd_len; ++i) {
		real r = x[XDIM][i] * x[XDIM][i];
		r += x[YDIM][i] * x[YDIM][i];
		r += x[ZDIM][i] * x[ZDIM][i];
		r = sqrt(r);
		u[D_i][i] = 1.0e-4;
		u[tau_i][i] = 1.0 * exp(-r / 0.05) + u[D_i][i];
		u[Sx_i][i] = u[Sy_i][i] = u[Sz_i][i] = 0.0;
	}
	return u;
}

euler_type::vector_type euler_type::set_outflow(const geo::direction& dir,
		vector_type&& u) {
	/*	for (integer d = 0; d != NDIM; ++d) {
	 if (dir[d] == -1) {
	 u[sxi + d] = min(u[sxi + d], simd_vector(0.0));
	 } else if (dir[d] == +1) {
	 u[sxi + d] = max(u[sxi + d], simd_vector(0.0));
	 }
	 }*/
	return std::move(u);
}

std::pair<euler_type::vector_type, simd_vector> euler_type::physical_flux(
		const vector_type& u, const vector_type& v, integer dim) {
	vector_type f(NF);
	simd_vector a;
	const auto& D = u[D_i];
	const auto& tau = u[tau_i];
	const auto& Sx = u[Sx_i];
	const auto& Sy = u[Sy_i];
	const auto& Sz = u[Sz_i];

	const auto& rho = v[rho_i];
	const auto& eps = v[eps_i];
	const auto& uel = v[ux_i + dim];
	const auto& ux = v[ux_i];
	const auto& uy = v[uy_i];
	const auto& uz = v[uz_i];

	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto Winv = one / W;
	const auto vel = uel * Winv;
	const auto v2 = u2 * Winv * Winv;
	const auto p = (fgamma - one) * rho * eps;
	const auto h = rho * (one + eps) + p;
	const auto hinv = one / h;
	const auto c2 = fgamma * p * h;
	const auto vel2 = vel * vel;
	f[D_i] = D * vel;
	f[Sx_i] = Sx * vel;
	f[Sy_i] = Sy * vel;
	f[Sz_i] = Sz * vel;
	f[Sx_i + dim] += p;
	f[tau_i] = (tau + p) * vel;
	const auto onemv2c2 = one - v2 * c2;
	const auto onemv2c2inv = one / onemv2c2;
	const auto onemc2 = one - c2;
	const auto a1 = onemv2c2inv * onemc2 * vel;
	const auto a2 = onemv2c2inv * sqrt(c2 * Winv * (onemv2c2 - onemc2 * vel2));
	a = abs(a1) + a2;
	std::pair<vector_type, simd_vector> pr;
	pr.first = std::move(f);
	pr.second = std::move(a);
	return pr;

}

euler_type::vector_type euler_type::to_con(const vector_type& v) {
	std::vector<simd_vector> u(NF);

	auto& D = u[D_i];
	auto& tau = u[tau_i];
	auto& Sx = u[Sx_i];
	auto& Sy = u[Sy_i];
	auto& Sz = u[Sz_i];

	const auto& rho = v[rho_i];
	const auto& eps = v[eps_i];
	const auto& ux = v[ux_i];
	const auto& uy = v[uy_i];
	const auto& uz = v[uz_i];

	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto p = (fgamma - 1.0) * rho * eps;
	const auto h = rho * (one + fgamma * eps);
	const auto hW = W * h;
	D = W * rho;
	tau = hW * W - p;
	Sx = hW * ux;
	Sy = hW * uy;
	Sz = hW * uz;

	return u;
}

euler_type::vector_type euler_type::to_prim(const vector_type& u) {
	std::vector<simd_vector> v(NF);
	static const simd_vector k((fgamma - one) / fgamma);

	const auto& D = u[D_i];
	auto tau = u[tau_i];
	const auto& Sx = u[Sx_i];
	const auto& Sy = u[Sy_i];
	const auto& Sz = u[Sz_i];

	auto& rho = v[rho_i];
	auto& eps = v[eps_i];
	auto& ux = v[ux_i];
	auto& uy = v[uy_i];
	auto& uz = v[uz_i];

	simd_vector W, vel, f, dfdv;
	simd_vector S(sqrt(Sx * Sx + Sy * Sy + Sz * Sz));

	tau = max(tau, sqrt(S * S + D * D));

	vel = S / D;
	const auto Sinv = S.inv_or_zero();
	vel = min(abs(vel), simd_vector(0.999));
	for (int i = 0; i < 8; ++i) {
		W = sqrt(one / (one - vel * vel));
		const simd_vector Winv(one / W);
		const simd_vector W2inv(Winv * Winv);
		f = (one - k * W2inv) * S + vel * k * D * Winv - vel * tau;
		dfdv = k * (two * S * vel + D * (one - two * vel * vel) * W) - tau;
		const auto v0 = vel;
		vel = vel - f / dfdv;
		const simd_vector b(0.1);
		vel = max(b * v0, vel);
		vel = min(one - b * (one - v0), vel);
	}
	W = sqrt(one / (one - vel * vel));
	const auto SinvvelW = Sinv * vel * W;
	rho = D / W;
	eps = (tau / rho - W * W) / (fgamma * W * W - (fgamma - one));
	ux = Sx * SinvvelW;
	uy = Sy * SinvvelW;
	uz = Sz * SinvvelW;
	eps = max(eps, zero);
	return v;
}

euler_type::vector_type euler_type::to_output(const vector_type& u) {
	auto v = to_prim(u);
	auto& ux = v[ux_i];
	auto& uy = v[uy_i];
	auto& uz = v[uz_i];
	const auto u2 = ux * ux + uy * uy + uz * uz;
	const auto W = sqrt(one + u2);
	const auto Winv = one / W;
	ux *= Winv;
	uy *= Winv;
	uz *= Winv;
	return v;
}

std::pair<euler_type::vector_type, simd_vector> euler_type::to_fluxes(
		const vector_type& vl, const vector_type& vr, integer dim) {
	std::vector<simd_vector> flux(NF);
	std::pair<euler_type::vector_type, simd_vector> r;

	const auto ur = to_con(vr);
	const auto ul = to_con(vl);

	const auto fr = physical_flux(ur, vr, dim);
	const auto fl = physical_flux(ul, vl, dim);

	const auto a = max(fr.second, fl.second);
	for (integer f = 0; f != NF; ++f) {
		flux[f] = half * (fr.first[f] + fl.first[f]);
		flux[f] -= half * a * (ur[f] - ul[f]);
	}
	r.first = std::move(flux);
	r.second = std::move(a);
	return r;
}

bool euler_type::refinement_test(integer level,
		const std::vector<simd_vector>& x, const vector_type& u,
		const std::array<vector_type, NDIM>& dudx) {
	const auto max_level = opts.max_level;
	return level < max_level;
}

euler_type::vector_type euler_type::explicit_source(const vector_type&) {
	std::vector<simd_vector> s(NF);
	std::fill(s.begin(), s.end(), simd_vector(0));
	return s;
}

euler_type::vector_type euler_type::implicit_source(const vector_type&, real) {
	std::vector<simd_vector> s(NF);
	std::fill(s.begin(), s.end(), simd_vector(0));
	return s;
}

std::vector<std::string> euler_type::field_names() {
	std::vector<std::string> names(NF);
	names[D_i] = "rho";
	names[tau_i] = "eps";
	names[Sx_i] = "vx";
	names[Sy_i] = "vy";
	names[Sz_i] = "vz";
	return names;
}

