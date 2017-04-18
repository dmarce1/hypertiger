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
class euler_type {
public:
	static constexpr integer NF = 5;
	using vector_type = std::vector<simd_vector>;
	static vector_type initial_value(const std::vector<simd_vector>&, real);
	static vector_type set_outflow(const geo::direction& dir, vector_type&& u);
	static bool refinement_test(integer, const std::vector<simd_vector>&,
			const vector_type&, const std::array<vector_type, NDIM>&);
	static vector_type to_prim(const vector_type&);
	static vector_type to_con(const vector_type&);
	static std::pair<vector_type, simd_vector> physical_flux( const vector_type& , const vector_type&, integer dim);
	static std::pair<vector_type, simd_vector> to_fluxes(const vector_type&,
			const vector_type&, integer dim);
	static vector_type explicit_source(const vector_type&);
	static vector_type implicit_source(const vector_type&, real);
	static std::vector<std::string> field_names();
};

using physics = euler_type;

#endif /* PROBLEM_HPP_ */
