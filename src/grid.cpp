#include "future.hpp"
#include "grid.hpp"
#include "problem.hpp"
#include "options.hpp"
#include "node_server.hpp"
#include <silo.h>

#include <array>
#include <cmath>
#include <cassert>

#include <hpx/include/runtime.hpp>
#include <hpx/lcos/broadcast.hpp>
#include <hpx/runtime/threads/run_as_os_thread.hpp>

std::array<simd_vector, NDIM> grid::outflow_mask;

struct tls_data_t {
	std::vector<std::vector<real>> v;
	std::vector<std::vector<std::vector<real>>>dvdx;
	std::vector<std::vector<std::vector<real>>> dudx;
	std::vector<std::vector<std::vector<real>>> uf;
	std::vector<std::vector<real>> zz;
};

void grid::set_hydro_boundary(const std::vector<simd_vector>& data,
		const geo::direction& dir) {

	std::array<integer, NDIM> lb, ub;
	get_boundary_size(lb, ub, dir, OUTER, INX / 2, 1);
	integer iter = 0;

	for (integer i = lb[XDIM]; i < ub[XDIM]; ++i) {
		for (integer j = lb[YDIM]; j < ub[YDIM]; ++j) {
			for (integer k = lb[ZDIM]; k < ub[ZDIM]; ++k) {
				const integer iii = icoarse(i, j, k);
				for (integer f = 0; f != physics::NF; ++f) {
					U(iii, f) = data[iter];
					++iter;
				}
			}
		}

	}
}

std::vector<simd_vector> grid::get_hydro_boundary(const geo::direction& dir) {

	std::array<integer, NDIM> lb, ub;
	std::vector<simd_vector> data;
	integer size;
	size = physics::NF * get_boundary_size(lb, ub, dir, INNER, INX / 2, 1);
	data.resize(size);
	integer iter = 0;

	for (integer i = lb[XDIM]; i < ub[XDIM]; ++i) {
		for (integer j = lb[YDIM]; j < ub[YDIM]; ++j) {
			for (integer k = lb[ZDIM]; k < ub[ZDIM]; ++k) {
				for (integer f = 0; f != physics::NF; ++f) {
					const integer iii = icoarse(i, j, k);
					data[iter] = U(iii, f);
					++iter;
				}
			}
		}

	}
	return data;

}

std::vector<real> grid::get_flux_restrict(const std::array<integer, NDIM>& lb,
		const std::array<integer, NDIM>& ub, const geo::dimension& dim) const {

	std::vector<real> data;
	integer size = 1;
	for (auto& dim : geo::dimension::full_set()) {
		size *= (ub[dim] - lb[dim]);
	}
	size /= (NCHILD / 2);
	size *= physics::NF;
	data.reserve(size);
	for (integer i = lb[XDIM]; i < ub[XDIM]; i += 2) {
		for (integer j = lb[YDIM]; j < ub[YDIM]; j += 2) {
			for (integer k = lb[ZDIM]; k < ub[ZDIM]; k += 2) {
				for (integer f = 0; f != physics::NF; ++f) {
					real value = F[dim](i, j, k, f);
					if (dim == XDIM) {
						value += F[dim](i, j + 1, k + 1, f);
					} else if (dim == YDIM) {
						value += F[dim](i + 1, j, k + 1, f);
					} else {
						value += F[dim](i + 1, j + 1, k, f);
					}
					if (dim != XDIM) {
						value += F[dim](i + 1, j, k, f);
					}
					if (dim != YDIM) {
						value += F[dim](i, j + 1, k, f);
					}
					if (dim != ZDIM) {
						value += F[dim](i, j, k + 1, f);
					}
					value /= real(4);
					data.push_back(value);
				}
			}
		}
	}
	return data;
}

void grid::set_flux_restrict(const std::vector<real>& data,
		const std::array<integer, NDIM>& lb,
		const std::array<integer, NDIM>& ub, const geo::dimension& dim) {

	integer index = 0;
	for (integer i = lb[XDIM]; i < ub[XDIM]; ++i) {
		for (integer j = lb[YDIM]; j < ub[YDIM]; ++j) {
			for (integer k = lb[ZDIM]; k < ub[ZDIM]; ++k) {
				for (integer f = 0; f != physics::NF; ++f) {
					F[dim](i, j, k, f) = data[index];
					++index;
				}
			}
		}
	}
}

void grid::set_prolong(const std::vector<simd_vector>& data) {

	integer index = 0;
	for (integer i = 1; i != NX / 2 - 1; ++i) {
		for (integer j = 1; j != NX / 2 - 1; ++j) {
			for (integer k = 1; k != NX / 2 - 1; ++k) {
				const integer iii = icoarse(i, j, k);
				for (integer field = 0; field != physics::NF; ++field) {
					U(iii, field) = data[index];
					++index;
				}
			}
		}
	}
}

std::vector<simd_vector> grid::get_prolong(const std::array<integer, NDIM>& lb,
		const std::array<integer, NDIM>& ub) {

	std::vector<simd_vector> data;

	integer size = physics::NF;
	for (integer dim = 0; dim != NDIM; ++dim) {
		size *= (ub[dim] - lb[dim]);
	}
	data.reserve(size / NCHILD);
	auto lb0 = lb;
	auto ub0 = ub;
	for (integer d = 0; d != NDIM; ++d) {
		lb0[d] /= 2;
		ub0[d] = (ub[d] - 1) / 2 + 1;
	}

	for (integer i = lb[XDIM]; i != ub[XDIM]; ++i) {
		const real xsgn = (i % 2) ? +1 : -1;
		for (integer j = lb[YDIM]; j != ub[YDIM]; ++j) {
			const real ysgn = (j % 2) ? +1 : -1;
#pragma GCC ivdep
			for (integer k = lb[ZDIM]; k != ub[ZDIM]; ++k) {
				for (integer f = 0; f != physics::NF; ++f) {
					data.push_back(simd_vector(U(i, j, k, f)));
				}
			}
		}
	}

	return data;
}

std::vector<real> grid::get_restrict() const {

	constexpr
	integer Size = physics::NF * INX * INX * INX / NCHILD;
	std::vector<real> data;
	data.reserve(Size);
	for (integer i = 1; i < NX / 2 - 1; ++i) {
		for (integer j = 1; j < NX / 2 - 1; ++j) {
			for (integer k = 1; k < NX / 2 - 1; ++k) {
				const integer iii = icoarse(i, j, k);
				for (integer field = 0; field != physics::NF; ++field) {
					data.push_back(U(iii, field).sum() / real(NCHILD));
				}
			}
		}
	}
	return data;
}

void grid::set_restrict(const std::vector<real>& data,
		const geo::octant& octant) {

	integer index = 0;
	for (integer i = BW; i != NX / 2; ++i) {
		for (integer j = BW; j != NX / 2; ++j) {
			for (integer k = BW; k != NX / 2; ++k) {
				for (integer field = 0; field != physics::NF; ++field) {
					const integer i0 = i + octant[XDIM] * INX / 2;
					const integer j0 = j + octant[YDIM] * INX / 2;
					const integer k0 = k + octant[ZDIM] * INX / 2;
					U(i0, j0, k0, field) = data[index];
					++index;
				}
			}
		}
	}
	assert(index == data.size());
}

bool grid::refine_me(integer lev) const {
	bool rc;
	if (lev < 1) {
		rc = true;
	} else if (lev >= opts.max_level) {
		rc = false;
	} else {
		std::vector<simd_vector> u(physics::NF);
		std::array<std::vector<simd_vector>, NDIM> dudx;
		for (integer d = 0; d != NDIM; ++d) {
			dudx[d].resize(physics::NF);
		}
		for (integer i = 1; i != NX / 2 - 1 && !rc; ++i) {
			for (integer j = 1; j != NX / 2 - 1 && !rc; ++j) {
				for (integer k = 1; k != NX / 2 - 1 && !rc; ++k) {
					const integer iii = icoarse(i, j, k);
					u = U[iii];
					for (integer d = 0; d != NDIM; ++d) {
						const auto up = U.get_shift_vec(d, +1, iii);
						const auto um = U.get_shift_vec(d, -1, iii);
						for (integer f = 0; f != physics::NF; ++f) {
							dudx[d][f] = (up[f] - um[f])
									* simd_vector(0.5 / dx);
						}
					}
					const auto test = physics::refinement_test;
					if (test(lev, X[iii], U[iii], dudx)) {
						rc = true;
						break;
					}
				}
			}
		}
	}
	return rc;
}

void grid::set_coordinates() {
	static std::once_flag flag;
	for (integer i = 0; i != NX; ++i) {
		for (integer j = 0; j != NX; ++j) {
			for (integer k = 0; k != NX; ++k) {
				X(i, j, k, XDIM) = (real(i - BW) + HALF) * dx + xmin[XDIM];
				X(i, j, k, YDIM) = (real(j - BW) + HALF) * dx + xmin[YDIM];
				X(i, j, k, ZDIM) = (real(k - BW) + HALF) * dx + xmin[ZDIM];
			}
		}
	}
	is_physical[FXM] = bool((xmin[XDIM] - dx) < 0.0);
	is_physical[FYM] = bool((xmin[YDIM] - dx) < 0.0);
	is_physical[FZM] = bool((xmin[ZDIM] - dx) < 0.0);
	is_physical[FXP] = bool((xmin[XDIM] + dx) > opts.xscale);
	is_physical[FYP] = bool((xmin[YDIM] + dx) > opts.xscale);
	is_physical[FZP] = bool((xmin[ZDIM] + dx) > opts.xscale);
	std::call_once(flag, [this]() {
		for (integer f = 0; f != NDIM; ++f) {
			outflow_mask[NDIM] = simd_vector(0.0);
		}
		for (integer i = 0; i != 2; ++i) {
			for (integer j = 0; j != 2; ++j) {
				outflow_mask[XDIM][geo::octant( {0, i, j})] = ONE;
				outflow_mask[YDIM][geo::octant( {i, 0, j})] = ONE;
				outflow_mask[ZDIM][geo::octant( {i, j, 0})] = ONE;
			}
		}
	});
}

void grid::allocate() {
	set_coordinates();
}

grid::grid() {
}

grid::grid(const real _dx, std::array<real, NDIM> _xmin, bool initialize) {
	dx = _dx;
	xmin = _xmin;
	allocate();
	if (initialize) {
		for (integer i = 0; i != NX / 2; ++i) {
			for (integer j = 0; j != NX / 2; ++j) {
				for (integer k = 0; k != NX / 2; ++k) {
					const integer iii = icoarse(i, j, k);
					const auto x = X[iii];
					auto u = physics::initial_value(x, dx);
					U.set(u, iii);
				}
			}
		}
	}
}

simd_grid<physics::NF> grid::primitives() const {
	simd_grid<physics::NF> V;
	for (integer i = 0; i != V_N3; ++i) {
		const auto uin = U[i];
		const auto vout = physics::to_prim(uin);
		V.set(vout, i);
	}
	return V;
}

simd_grid<physics::NF> grid::output_vars() const {
	simd_grid<physics::NF> V;
	for (integer i = 0; i != V_N3; ++i) {
		const auto uin = U[i];
		const auto vout = physics::to_output(uin);
		V.set(vout, i);
	}
	return V;
}

real grid::compute_fluxes() {
	simd_vector max_lambda(0.0);
	simd_grid<physics::NF> VR;
	simd_grid<physics::NF> VL;
	const auto V = primitives();
	for (integer dim = 0; dim != NDIM; ++dim) {
		for (integer i = 0; i != U.size(); ++i) {
			const auto sp = V.get_shift(dim, +1, i) - V(i);
			const auto sm = V(i) - V.get_shift(dim, -1, i);
			auto s0 = min(abs(sp), abs(sm));
			s0 = (copysign(s0, sp) + copysign(s0, sm)) * 0.25;
			VR(i) = V(i) + s0;
			VL(i) = V(i) - s0;
		}
		for (integer i = 1; i != NX / 2; ++i) {
			for (integer j = 1; j != NX / 2; ++j) {
				for (integer k = 1; k != NX / 2; ++k) {
					const integer iii = icoarse(i, j, k);
					const auto left = VR.get_shift_vec(dim, -1, iii);
					const auto right = VL[iii];
					const auto f = physics::to_fluxes(left, right, dim);
					F[dim].set(std::move(f.first), iii);
					max_lambda = max(max_lambda, f.second);
				}
			}
		}
	}
	return max_lambda.max();
}

void grid::store() {
	U0 = U;
}

void grid::set_physical_boundaries(const geo::face& face, real t) {
	const auto dim = face.get_dimension();
	const auto side = face.get_side();
	const integer dnk = DN[dim];
	const integer klb = side == geo::MINUS ? 0 : NX - BW;
	const integer kub = side == geo::MINUS ? BW : NX;
	const integer ilb = 0;
	const integer iub = NX;
	const integer jlb = 0;
	const integer jub = NX;

	integer i, j, k, k0;
	const integer& x = dim == XDIM ? k : i;
	const integer& y = dim == YDIM ? k : (dim == XDIM ? i : j);
	const integer& z = dim == ZDIM ? k : j;
	const integer& x0 = dim == XDIM ? k0 : x;
	const integer& y0 = dim == YDIM ? k0 : y;
	const integer& z0 = dim == ZDIM ? k0 : z;

	for (integer field = 0; field != physics::NF; ++field) {
		for (k = klb; k != kub; ++k) {
			for (j = jlb; j != jub; ++j) {
				for (i = ilb; i != iub; ++i) {
					k0 = side == geo::MINUS ? BW : NX - BW - 1;
					U(x, y, z, field) = U(x0, y0, z0, field);
				}
			}
		}
	}
}

void grid::compute_sources(real t) {
	for (integer i = 0; i != V_N3; ++i) {
		const auto u = U[i];
		const auto s = physics::explicit_source(u);
		dUdt.set(s, i);
	}
}

void grid::compute_dudt() {
	const auto dxinv = simd_vector(1.0) / simd_vector(dx);
	oflux.resize(physics::NF);
	std::vector<simd_vector> outflow(physics::NF);
	for (integer dim = 0; dim != NDIM; ++dim) {
		for (integer i = 1; i != NX / 2 - 1; ++i) {
			for (integer j = 1; j != NX / 2 - 1; ++j) {
				for (integer k = 1; k != NX / 2 - 1; ++k) {
					for (integer f = 0; f != physics::NF; ++f) {
						const integer iii = f + physics::NF * icoarse(i, j, k);
						auto tmp = F[dim](iii);
						tmp -= F[dim].get_shift(dim, +1, iii);
						tmp *= simd_vector(dxinv);
						dUdt(iii) += tmp;
					}
				}
			}
		}
	}
	std::fill(oflux.begin(), oflux.end(), simd_vector(0.0));
	for (integer j = 1; j != NX / 2 - 1; ++j) {
		for (integer k = 1; k != NX / 2 - 1; ++k) {
			for (integer f = 0; f != physics::NF; ++f) {
				if (is_physical[FXM]) {
					const integer i = f + physics::NF * icoarse(1, j, k);
					oflux[f] += F[XDIM](i) * outflow_mask[XDIM];
				}
				if (is_physical[FXP]) {
					const integer i = physics::NF * icoarse(NX / 2 - 1, j, k);
					oflux[f] += F[XDIM](i + f) * outflow_mask[XDIM];
				}
				if (is_physical[FYM]) {
					const integer i = f + physics::NF * icoarse(j, 1, k);
					oflux[f] += F[YDIM](i) * outflow_mask[YDIM];
				}
				if (is_physical[FYP]) {
					const integer i = physics::NF * icoarse(j, NX / 2 - 1, k);
					oflux[f] += F[YDIM](i + f) * outflow_mask[YDIM];
				}
				if (is_physical[FZM]) {
					const integer i = f + physics::NF * icoarse(j, k, 1);
					oflux[f] += F[ZDIM](i) * outflow_mask[ZDIM];
				}
				if (is_physical[FZP]) {
					const integer i = physics::NF * icoarse(j, k, NX / 2 - 1);
					oflux[f] += F[ZDIM](i + f) * outflow_mask[ZDIM];
				}
			}
		}
	}
}

std::vector<simd_vector> grid::next_u(integer rk, real t, real dt,
		std::vector<simd_vector>&& O) {
	const simd_vector beta(rk == 0 ? 1.0 : 0.5);
	const simd_vector onembeta = simd_vector(1.0) - beta;
	for (integer i = 0; i != U.size(); ++i) {
		const auto U1 = U(i) + dUdt(i) * simd_vector(dt);
		U(i) = onembeta * U0(i) + beta * U1;
	}
	const simd_vector dx2 = simd_vector(dx * dx);
	for (integer f = 0; f != physics::NF; ++f) {
		O[f] = beta * (O[f] + oflux[f] * dx2);
	}
	return O;
}

inline bool float_eq(xpoint_type a, xpoint_type b) {
	constexpr static xpoint_type eps = 0.00000011920928955078125; // std::pow(xpoint_type(2), -23);
// 	const xpoint_type eps = std::pow(xpoint_type(2), -23);
	return std::abs(a - b) < eps;
}

bool grid::xpoint_eq(const xpoint& a, const xpoint& b) {
	bool rc = true;
	for (integer d = 0; d != NDIM; ++d) {
		if (!float_eq(a[d], b[d])) {
			rc = false;
			break;
		}
	}
	return rc;
}

bool grid::node_point::operator==(const node_point& other) const {
	return xpoint_eq(other.pt, pt);
}

bool grid::node_point::operator<(const node_point& other) const {
	bool rc = false;
	for (integer d = 0; d != NDIM; ++d) {
		if (!float_eq(pt[d], other.pt[d])) {
			rc = (pt[d] < other.pt[d]);
			break;
		}
	}
	return rc;
}

void grid::merge_output_lists(grid::output_list_type& l1,
		grid::output_list_type&& l2) {

	std::unordered_map<zone_int_type, zone_int_type> index_map;

	if (l2.zones.size() > l1.zones.size()) {
		std::swap(l1, l2);
	}
	for (auto i = l2.nodes.begin(); i != l2.nodes.end(); ++i) {
		zone_int_type index, oindex;
		auto this_x = *i;
		oindex = this_x.index;
		auto j = l1.nodes.find(this_x);
		if (j != l1.nodes.end()) {
			index = j->index;
		} else {
			index = l1.nodes.size();
			this_x.index = index;
			l1.nodes.insert(this_x);
		}
		index_map[oindex] = index;
	}
	integer zzz = l1.zones.size();
	l1.zones.resize(zzz + l2.zones.size());
	for (auto i = l2.zones.begin(); i != l2.zones.end(); ++i) {
		l1.zones[zzz] = index_map[*i];
		++zzz;
	}
	for (integer field = 0; field < physics::NF; ++field) {
		const auto l1sz = l1.data[field].size();
		l1.data[field].resize(l1sz + l2.data[field].size());
		std::move(l2.data[field].begin(), l2.data[field].end(),
				l1.data[field].begin() + l1sz);
	}
}

grid::output_list_type grid::get_output_list() const {
	output_list_type rc;
	const integer vertex_order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };

	std::set<node_point>& node_list = rc.nodes;
	std::vector<zone_int_type>& zone_list = rc.zones;
	std::array<std::vector<real>, physics::NF> &data = rc.data;

	const auto V = output_vars();

	zone_list.reserve(cube(NX - 2 * BW) * NCHILD);
	for (integer i = BW; i != NX - BW; ++i) {
		for (integer j = BW; j != NX - BW; ++j) {
			for (integer k = BW; k != NX - BW; ++k) {
				for (integer ci = 0; ci != NVERTEX; ++ci) {
					const integer vi = vertex_order[ci];
					const integer xi = (vi >> 0) & 1;
					const integer yi = (vi >> 1) & 1;

					const integer zi = (vi >> 2) & 1;
					node_point this_x;
					this_x.pt[XDIM] = X(i, j, k, XDIM) + (real(xi) - HALF) * dx;
					this_x.pt[YDIM] = X(i, j, k, YDIM) + (real(yi) - HALF) * dx;
					this_x.pt[ZDIM] = X(i, j, k, ZDIM) + (real(zi) - HALF) * dx;
					auto iter = node_list.find(this_x);
					integer index;
					if (iter != node_list.end()) {
						index = iter->index;
					} else {
						index = node_list.size();
						this_x.index = index;
						node_list.insert(this_x);
					}
					zone_list.push_back(index);
				}
				for (integer field = 0; field != physics::NF; ++field) {
					data[field].push_back(V(i, j, k, field));
				}
			}
		}
	}

	return rc;
}

void make_names(std::vector<char*>& names, std::vector<int>& types,
		std::string dirname, std::string base, std::string title, int type) {
	const integer sz = names.size();
	for (integer i = 0; i != sz; ++i) {
		std::string name = dirname + base + std::string(".") + std::to_string(i)
				+ std::string(".silo:") + title;
		names[i] = (char*) malloc(name.size() + 1);
		types[i] = type;
		strcpy(names[i], name.c_str());
	}
}

void delete_names(std::vector<char*>& names) {
	const integer sz = names.size();
	for (integer i = 0; i != sz; ++i) {
		free(names[i]);
	}
}

void grid::output_header(std::string dirname, std::string base, real t,
		int cycle, int procs) {
	hpx::threads::run_as_os_thread(
			[&]() {
				auto olist = DBMakeOptlist(1);
				double time = double(t);
				int ndim = 3;
				DBAddOption(olist, DBOPT_DTIME, &time);
				std::string filename = dirname + base + std::string(".silo");
				DBfile *db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
				assert(db);
				std::vector<char*> names(procs);
				std::vector<int> types(procs);
				make_names(names, types, dirname, base, "mesh", DB_UCDMESH);
				DBPutMultimesh(db, "mesh", procs, names.data(), types.data(), olist);
				delete_names(names);
				static const auto field_names = physics::field_names();
				for (int field = 0; field != physics::NF; ++field) {
					make_names(names, types, dirname, base, field_names[field], DB_UCDVAR);
					DBPutMultivar(db, field_names[field].c_str(), procs, names.data(), types.data(), olist);
					delete_names(names);
				}
				DBFreeOptlist(olist);
				DBClose(db);
			}).get();
}

void grid::output(const output_list_type& olists, std::string _dirname,
		std::string _base, real _t, int cycle) {
	hpx::threads::run_as_os_thread(
			[&](const std::string& dirname, const std::string& base, real t) {
				const std::set<node_point>& node_list = olists.nodes;
				const std::vector<zone_int_type>& zone_list = olists.zones;

				const int nzones = zone_list.size() / NVERTEX;
				std::vector<int> zone_nodes;
				zone_nodes = std::move(zone_list);

				const int nnodes = node_list.size();
				std::vector<double> x_coord(nnodes);
				std::vector<double> y_coord(nnodes);
				std::vector<double> z_coord(nnodes);
				std::array<double*, NDIM> node_coords = {x_coord.data(), y_coord.data(), z_coord.data()};
				for (auto iter = node_list.begin(); iter != node_list.end(); ++iter) {
					const integer i = iter->index;
					x_coord[i] = iter->pt[0];
					y_coord[i] = iter->pt[1];
					z_coord[i] = iter->pt[2];
				}

				constexpr int nshapes = 1;
				int shapesize[1] = {NVERTEX};
				int shapetype[1] = {DB_ZONETYPE_HEX};
				int shapecnt[1] = {nzones};
				const char* coord_names[NDIM] = {"x", "y", "z"};
				auto olist = DBMakeOptlist(1);
				double time = double(t);
				int ndim = 3;
				DBAddOption(olist, DBOPT_DTIME, &time);
				std::string filename = dirname + base;
				DBfile *db = DBCreateReal(filename.c_str(), DB_CLOBBER, DB_LOCAL, "Euler Mesh", DB_PDB);
				assert(db);
				DBPutZonelist2(db, "zones", nzones, int(NDIM), zone_nodes.data(), nzones * NVERTEX, 0, 0, 0, shapetype, shapesize,
						shapecnt, nshapes, olist);
				DBPutUcdmesh(db, "mesh", int(NDIM), const_cast<char**>(coord_names), node_coords.data(), nnodes, nzones, "zones", nullptr, DB_DOUBLE,
						olist);
				DBFreeOptlist(olist);
				for (int field = 0; field != physics::NF; ++field) {
					auto olist = DBMakeOptlist(1);
					double time = double(t);
					int istrue = 1;
					int isfalse = 0;
					DBAddOption(olist, DBOPT_DTIME, &time);
					static const auto field_names = physics::field_names();
					DBPutUcdvar1(db, field_names[field].c_str(), "mesh", const_cast<void*>(reinterpret_cast<const void*>(olists.data[field].data())), nzones, nullptr, 0, DB_DOUBLE, DB_ZONECENT,
							olist);
					DBFreeOptlist(olist);
				}
				DBClose(db);
			}, _dirname, _base, _t).get();
}

std::size_t grid::load(FILE* fp) {
	static hpx::mutex mtx;
	std::size_t cnt = 0;
	{
		static std::atomic<bool> statics_loaded(false);
		bool expected = false;
		if (statics_loaded.compare_exchange_strong(expected, true)) {
			cnt += std::fread(&opts.xscale, sizeof(real), 1, fp) * sizeof(real);
			statics_loaded = true;
		} else {
			std::size_t offset = sizeof(real);
			std::fseek(fp, offset, SEEK_CUR);
			cnt += offset;
		}
	}
	allocate();

	for (integer i = 1; i < NX / 2 - 1; ++i) {
		for (integer j = 1; j < NX / 2 - 1; ++j) {
			for (integer k = 1; k < NX / 2 - 1; ++k) {
				const integer ic = icoarse(i, j, k);
				for (integer f = 0; f != physics::NF; ++f) {
					auto& tmp = U(ic, f);
					cnt += std::fread(&tmp, sizeof(simd_vector), 1, fp)
							* sizeof(simd_vector);
				}
			}
		}
	}
	set_coordinates();
	return cnt;
}

std::size_t grid::save(std::ostream& strm) const {
	std::size_t cnt = 0;

	cnt += write(strm, opts.xscale);

	for (integer i = 1; i < NX / 2 - 1; ++i) {
		for (integer j = 1; j < NX / 2 - 1; ++j) {
			for (integer k = 1; k < NX / 2 - 1; ++k) {
				const integer ic = icoarse(i, j, k);
				for (integer f = 0; f != physics::NF; ++f) {
					auto tmp = U(ic, f);
					cnt += write(strm, &tmp, 1);
				}
			}
		}
	}
	return cnt;
}
