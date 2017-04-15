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
#include <hpx/include/lcos.hpp>

euler_type::vector_type euler_type::initial_value(
		const std::vector<simd_vector>& x, real) {
	std::vector<simd_vector> u(NF);
	for (integer i = 0; i != simd_len; ++i) {
		if (x[ZDIM][i] > 0.0) {
			u[d_i][i] = 1.0;
			u[e_i][i] = 1.0;
		} else {
			u[d_i][i] = 1.0e-1;
			u[e_i][i] = 1.25e-1;
		}
		u[sxi][i] = u[syi][i] = u[szi][i] = 0.0;
	}
	return u;
}

euler_type::vector_type euler_type::set_outflow(const geo::direction& dir,
		vector_type&& u) {
	for (integer d = 0; d != NDIM; ++d) {
		if (dir[d] == -1) {
			u[sxi + d] = min(u[sxi + d], simd_vector(0.0));
		} else if (dir[d] == +1) {
			u[sxi + d] = max(u[sxi + d], simd_vector(0.0));
		}
	}
	return std::move(u);
}

euler_type::vector_type euler_type::to_prim(const vector_type& u) {
	std::vector<simd_vector> v(NF);
	v[d_i] = u[d_i];

	/*********OPTIMIZE******/
	auto dinv = u[d_i];
	for (integer i = 0; i != simd_len; ++i) {
		if (dinv[i] != 0.0) {
			dinv[i] = 1.0 / dinv[i];
		}
	}
	/***************/

	v[vxi] = u[sxi] * dinv;
	v[vyi] = u[syi] * dinv;
	v[vzi] = u[szi] * dinv;
	v[t_i] = u[e_i];
	v[t_i] -= simd_vector(0.5)
			* (u[vxi] * u[vxi] + u[vyi] * u[vyi] + u[vzi] * u[vzi]);
	v[t_i] = max(simd_vector(fgamma - 1.0) * v[t_i] * dinv, simd_vector(0.0));
	return v;
}

std::pair<euler_type::vector_type, simd_vector> euler_type::to_fluxes(
		const vector_type& vl, const vector_type& vr, integer dim) {
	std::vector<simd_vector> f(NF);
	std::pair<euler_type::vector_type, simd_vector> r;

	/*********OPTIMIZE******/
	auto dr_inv = vr[d_i];
	auto dl_inv = vl[d_i];
	for (integer i = 0; i != simd_len; ++i) {
		if (dr_inv[i] != 0.0) {
			dr_inv[i] = 1.0 / dr_inv[i];
		}
		if (dl_inv[i] != 0.0) {
			dl_inv[i] = 1.0 / dl_inv[i];
		}
	}
	/***************/

	f[d_i] = simd_vector(0.5)
			* (vl[d_i] * vl[vxi + dim] + vr[d_i] * vr[vxi + dim]);
	for (integer d2 = 0; d2 != NDIM; ++d2) {
		f[sxi + d2] = +vr[d_i] * vr[vxi + dim] * vr[vxi + d2];
		f[sxi + d2] += vl[d_i] * vl[vxi + dim] * vl[vxi + d2];
	}
	const auto pl = vl[d_i] * vl[t_i];
	const auto pr = vr[d_i] * vr[t_i];
	const auto hl = pl * simd_vector(fgamma / (fgamma - 1.0));
	const auto hr = pr * simd_vector(fgamma / (fgamma - 1.0));
	const auto cl = sqrt(simd_vector(fgamma) * pl * dl_inv);
	const auto cr = sqrt(simd_vector(fgamma) * pr * dr_inv);
	const auto al = abs(vl[vxi + dim]) + cl;
	const auto ar = abs(vr[vxi + dim]) + cr;

	f[sxi + dim] += pl + pr;
	f[e_i] = simd_vector(0.5) * (vl[vxi + dim] * hl + vr[vxi + dim] * hr);
	for (integer d2 = 0; d2 != NDIM; ++d2) {
		f[sxi + d2] *= simd_vector(0.5);
	}
	auto a = max(al, ar);
	f[d_i] -= simd_vector(0.5) * a * (vr[d_i] - vl[d_i]);
	f[sxi] -= simd_vector(0.5) * a * (vr[d_i] * vr[vxi] - vl[d_i] * vl[vxi]);
	f[syi] -= simd_vector(0.5) * a * (vr[d_i] * vr[vyi] - vl[d_i] * vl[vyi]);
	f[szi] -= simd_vector(0.5) * a * (vr[d_i] * vr[vzi] - vl[d_i] * vl[vzi]);
	f[e_i] -= simd_vector(0.5) * a * (vr[d_i] * vr[t_i] - vl[d_i] * vl[t_i])
			/ simd_vector(fgamma - 1.0);
	r.first = std::move(f);
	r.second = std::move(a);
	return r;
}

bool euler_type::refinement_test(integer level,
		const std::vector<simd_vector>& x, const vector_type& u,
		const std::array<vector_type, NDIM>& dudx) {
	const auto max_level = opts.xscale;
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
	names[d_i] = "d";
	names[e_i] = "e";
	names[sxi] = "sx";
	names[syi] = "sy";
	names[szi] = "sz";
	return names;
}

