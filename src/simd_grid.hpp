/*
 * simd_grid.hpp
 *
 *  Created on: Apr 12, 2017
 *      Author: dminacore
 */

#ifndef SRC_SIMD_GRID_HPP_
#define SRC_SIMD_GRID_HPP_

#include "defs.hpp"
#include <cassert>
#include "simd.hpp"
#include "vector"
#include <cstring>

constexpr integer V_N3 = N3 / simd_len;

template<integer NF>
class simd_grid {
private:
	std::vector<simd_vector> data;
	using iterator_type = typename std::vector<simd_vector>::iterator;
	static integer index(integer i, integer f) {
		return i * NF + f;
	}
public:
	inline constexpr std::size_t size() const {
		return NF * V_N3;
	}
	inline simd_grid() :
			data(size()) {
	}
	inline iterator_type begin() {
		return data.begin();
	}
	inline iterator_type end() {
		return data.end();
	}
	inline real operator()(integer i, integer j, integer k, integer f) const {
		const integer i1 = (k / 2) * (NX * NX / 4) + (j / 2) * (NX / 2) + i / 2;
		const integer i2 = (k % 2) * 4 + (j % 2) * 2 + i % 2;
		return data[i1 * NF + f][i2];
	}
	inline real& operator()(integer i, integer j, integer k, integer f) {
		const integer i1 = (k / 2) * (NX * NX / 4) + (j / 2) * (NX / 2) + i / 2;
		const integer i2 = (k % 2) * 4 + (j % 2) * 2 + i % 2;
		return data[i1 * NF + f][i2];
	}
	inline std::vector<simd_vector> operator[](integer i) const {
		std::vector<simd_vector> v(NF);
		const auto* base = data.data() + i * NF;
		const auto len = sizeof(simd_vector) * NF;
		std::memcpy(v.data(), base, len);
		return v;
	}
	inline void set(const std::vector<simd_vector>& v, integer i) {
		auto* base = data.data() + i * NF;
		const auto len = sizeof(simd_vector) * NF;
		std::memcpy(base, v.data(), len);
	}
	inline simd_vector operator()(integer i, integer f) const {
		return data[index(i, f)];
	}
	inline simd_vector operator()(integer i) const {
		return data[i];
	}
	inline simd_vector& operator()(integer i) {
		return data[i];
	}
	inline simd_vector& operator()(integer i, integer f) {
		return data[index(i, f)];
	}
	inline simd_vector get_shift(const integer dim, const integer sign,
			const integer i1) const {
		constexpr integer DN2[NDIM] = { 1, NX / 2, NX * NX / 4 };
		const integer shift = NF * std::copysign(DN2[dim], sign);
		const integer i2 = (i1 + shift + size()) % size();
		simd_vector a = (*this)(i1).shuffle(dim);
		simd_vector b = (*this)(i2).shuffle(dim);
		if (sign > 0) {
			return blend(a, b, dim);
		} else {
			return blend(b, a, dim);
		}
	}
	inline std::vector<simd_vector> get_shift_vec(integer dim, integer sign,
			integer i) const {
		std::vector<simd_vector> v(NF);
		for (integer f = 0; f != NF; ++f) {
			v[f] = get_shift(dim, sign, i * NF + f);
		}
		return v;
	}
	template<class Arc>
	inline void serialize(Arc& arc, const unsigned) {
		arc & data;
	}
};

#endif /* SRC_SIMD_GRID_HPP_ */
