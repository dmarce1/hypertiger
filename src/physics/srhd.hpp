/*
 * problem.hpp
 *
 *  Created on: May 29, 2015
 *      Author: dmarce1
 */

#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_

#include "defs.hpp"
#include "simd.hpp"
#include <vector>
#include <functional>
#include "geometry.hpp"

class srhd {
public:
	static constexpr integer NF = 5;
	using vector_type = std::array<simd_vector,NF>;
	static void initial_value(vector_type&,
			const std::array<simd_vector, NDIM>&, real);
	static void set_outflow(vector_type&, const geo::direction& dir);
	static bool refinement_test(integer, const std::array<simd_vector, NDIM>&,
			const vector_type&, const std::array<vector_type, NDIM>&);
	static void to_output(vector_type&, const vector_type&);
	static void to_prim(vector_type&, const vector_type&);
	static void to_con(vector_type&, const vector_type&);
	static void physical_flux(vector_type&, simd_vector&, const vector_type&,
			const vector_type&, integer dim);
	static void to_fluxes(vector_type&, simd_vector&, const vector_type&,
			const vector_type&, integer dim);
	static void explicit_source(vector_type& s, const vector_type& u, const vector_type& v);
	static void implicit_source(vector_type&, const vector_type&, real);
	static std::vector<std::string> field_names();
};

#endif /* PROBLEM_HPP_ */
