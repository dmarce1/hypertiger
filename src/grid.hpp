/*
 * grid.hpp
 *
 *  Created on: May 26, 2015
 *      Author: dmarce1
 */

#ifndef GRID_HPP_
#define GRID_HPP_

#include "simd_grid.hpp"
#include "simd.hpp"
#include "defs.hpp"
#include "space_vector.hpp"
#include "geometry.hpp"
#include "problem.hpp"
#include "real.hpp"

#include <hpx/runtime/serialization/serialize.hpp>
#include <hpx/runtime/serialization/set.hpp>
#include <hpx/runtime/serialization/array.hpp>
#include <hpx/runtime/serialization/vector.hpp>
#include <hpx/traits/is_bitwise_serializable.hpp>
#include <valarray>

#include <iostream>
#include <valarray>
#include <functional>
#include <list>
#include <memory>
#include <set>

class struct_eos;

typedef real xpoint_type;
typedef int zone_int_type;

class grid {
public:
	typedef std::array<xpoint_type, NDIM> xpoint;
	struct node_point;
private:
	static std::array<simd_vector, NDIM> outflow_mask;

	simd_grid<physics::NF> U;
	simd_grid<physics::NF> U0;
	simd_grid<physics::NF> dUdt;
	std::array<simd_grid<physics::NF>, NDIM> F;
	simd_grid<NDIM> X;
	real dx;
	std::array<real, NDIM> xmin;
	std::array<bool, NFACE> is_physical;
	std::vector<simd_vector> oflux;
	static bool xpoint_eq(const xpoint& a, const xpoint& b);
public:
	static void set_max_level(integer l);
	real get_dx() const {
		return dx;
	}
	simd_grid<physics::NF> primitives() const;
	void set_coordinates();
	void set_hydro_boundary(const std::vector<simd_vector>&,
			const geo::direction&);
	std::vector<simd_vector> get_hydro_boundary(const geo::direction& face);
	struct output_list_type;
	static void merge_output_lists(output_list_type& l1, output_list_type&& l2);
	std::vector<real> get_restrict() const;
	std::vector<real> get_flux_restrict(const std::array<integer, NDIM>& lb,
			const std::array<integer, NDIM>& ub, const geo::dimension&) const;
	std::vector<simd_vector> get_prolong(const std::array<integer, NDIM>& lb,
			const std::array<integer, NDIM>& ub);
	void set_prolong(const std::vector<simd_vector>&);
	void set_restrict(const std::vector<real>&, const geo::octant&);
	void set_flux_restrict(const std::vector<real>&,
			const std::array<integer, NDIM>& lb,
			const std::array<integer, NDIM>& ub, const geo::dimension&);
	bool refine_me(integer lev) const;
	void compute_dudt();
	integer get_step() const;
	grid(real dx, std::array<real, NDIM> xmin, bool = false);
	grid();
	~grid() {
	}
	grid(const grid&) = delete;
	grid(grid&&) = default;
	grid& operator=(const grid&) = delete;
	grid& operator=(grid&&) = default;
	std::pair<space_vector, space_vector> find_axis() const;
	space_vector get_cell_center(integer i, integer j, integer k);
	void allocate();
	void store();
	void restore();
	real compute_fluxes();
	void compute_sources(real t);
	void set_physical_boundaries(const geo::face&, real t);
	std::vector<simd_vector> next_u(integer rk, real t, real dt,
			std::vector<simd_vector>&&);
	static void output(const output_list_type&, std::string, std::string,
			real t, int cycle);
	static void output_header(std::string, std::string, real t, int cycle,
			int procs);
	template<class Archive>
	void serialize(Archive& arc, const unsigned) {
		arc & dx;
		arc & xmin;
		arc & U;
		allocate();
	}
	std::size_t load(FILE* fp);
	std::size_t save(std::ostream& strm) const;
	output_list_type get_output_list() const;

	friend class node_server;
};

struct grid::node_point {
	xpoint pt;
	integer index;
	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & pt;
		arc & index;
	}
	bool operator==(const grid::node_point& other) const;
	bool operator<(const grid::node_point& other) const;
};

namespace hpx {
namespace traits {
template<>
struct is_bitwise_serializable<grid::node_point> : std::true_type {
};
}
}

struct grid::output_list_type {
	std::set<node_point> nodes;
	std::vector<zone_int_type> zones;
	std::array<std::vector<real>, physics::NF> data;
	template<class Arc>
	void serialize(Arc& arc, unsigned int) {
		arc & nodes;
		arc & zones;
		arc & data;
	}
}
;

#endif /* GRID_HPP_ */
