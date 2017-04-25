/*
 * problem.hpp
 *
 *  Created on: May 29, 2015
 *      Author: dmarce1
 */

#ifndef GRHD_HPP_
#define GRHD_HPP_

#include "defs.hpp"
#include "simd.hpp"
#include <vector>
#include <functional>
#include "geometry.hpp"
#include "srhd.hpp"

namespace gr {
union z4_t;
}

class grhd {
public:
	static constexpr integer NF = 38 + srhd::NF;
	using vector_type = std::array<simd_vector,NF>;
private:
	static gr::z4_t& get_z4(vector_type& ref);
	static const gr::z4_t& get_z4(const vector_type& ref);
public:
	static void initial_value(vector_type&,
			const std::array<simd_vector, NDIM>&, real);
	static void set_outflow(vector_type&, const geo::face& dir, real dx);
	static bool refinement_test(integer, const std::array<simd_vector, NDIM>&,
			const vector_type&, const std::array<vector_type, NDIM>&);
	static void to_output(vector_type&, const vector_type&);
	static void to_prim(vector_type&, const vector_type&);
	static void to_con(vector_type&, const vector_type&);
	static void physical_flux(vector_type&, simd_vector&, const vector_type&,
			const vector_type&, integer dim,
			const std::array<simd_vector, NDIM>&, real t);
	static void to_fluxes(vector_type&, simd_vector&, vector_type&,
			vector_type&, integer dim, const std::array<simd_vector, NDIM>&, real t);
	static void explicit_source(vector_type& s, const vector_type& u,
			const vector_type& v, const std::array<simd_vector, NDIM>&, real t);
	static void implicit_source(vector_type&, const vector_type&, real);
	static std::vector<std::string> field_names();
};

#endif /* PROBLEM_HPP_ */
