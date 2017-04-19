/*
 * simd_vector.hpp
 *
 *  Created on: Jun 7, 2015
 *      Author: dmarce1
 */

#ifndef SIMD_VECTOR_HPP_
#define SIMD_VECTOR_HPP_
#include "defs.hpp"
#include <cstdlib>
#include <cstdio>
#include "immintrin.h"

constexpr std::size_t simd_len = 8;

#if !defined(HPX_HAVE_DATAPAR)

#ifdef USE_SIMD
#if !defined(__MIC__) && !defined(__AVX512F__)
#define SIMD_SIZE 2
#define __mxd __m256d
#define _mmx_set_pd(d)    _mm256_set_pd((d),(d),(d),(d))
#define _mmx_add_pd(a,b)  _mm256_add_pd((a),(b))
#define _mmx_sub_pd(a,b)  _mm256_sub_pd((a),(b))
#define _mmx_mul_pd(a,b)  _mm256_mul_pd((a),(b))
#define _mmx_div_pd(a,b)  _mm256_div_pd((a),(b))
#define _mmx_sqrt_pd(a)   _mm256_sqrt_pd(a)
#define _mmx_max_pd(a, b) _mm256_max_pd((a),(b))
#define _mmx_min_pd(a, b) _mm256_min_pd((a),(b))
#define _mmx_or_pd(a, b)  _mm256_or_pd((a),(b))
#define _mmx_and_pd(a, b) _mm256_and_pd((a),(b))
#else
//#warning "Compiling for Intel MIC"
#define SIMD_SIZE 1
#define __mxd __m512d
#define _mmx_set_pd(d)    _mm512_set_pd((d),(d),(d),(d),(d),(d),(d),(d))
#define _mmx_add_pd(a,b)  _mm512_add_pd((a),(b))
#define _mmx_sub_pd(a,b)  _mm512_sub_pd((a),(b))
#define _mmx_mul_pd(a,b)  _mm512_mul_pd((a),(b))
#define _mmx_div_pd(a,b)  _mm512_div_pd((a),(b))
#define _mmx_sqrt_pd(a)   _mm512_sqrt_pd(a)
#define _mmx_max_pd(a, b) _mm512_max_pd((a),(b))
#endif
#else
#define SIMD_SIZE 8
#define __mxd double
#define _mmx_set_pd(d)    (d)
#define _mmx_add_pd(a,b)  ((a)+(b))
#define _mmx_sub_pd(a,b)  ((a)-(b))
#define _mmx_mul_pd(a,b)  ((a)*(b))
#define _mmx_div_pd(a,b)  ((a)/(b))
#define _mmx_sqrt_pd(a)   sqrt(a)
#define _mmx_max_pd(a, b) std::max((a),(b))
#endif

class simd_vector {
private:
	std::array<__mxd, SIMD_SIZE> v;
public:
	inline simd_vector shuffle(integer dim) const {
		simd_vector r;
		switch (dim) {
		case XDIM:
			r.v[0] = _mm256_permute_pd(v[0], 0x5);
			r.v[1] = _mm256_permute_pd(v[1], 0x5);
			break;
		case YDIM:
			r.v[0] = _mm256_permute4x64_pd(v[0], 0x4E);
			r.v[1] = _mm256_permute4x64_pd(v[1], 0x4E);
			break;
		case ZDIM:
			r.v[0] = v[1];
			r.v[1] = v[0];
			break;
		}
		return r;
	}
	template<class Arc>
	void serialize(Arc& arc, const unsigned) {
		for (integer i = 0; i != SIMD_SIZE; ++i) {
			for (integer j = 0; j != simd_len; ++j) {
				arc & v[i][j];
			}
		}
	}

	inline simd_vector inv_or_zero() const {
		simd_vector c;
		__mmask8 b;
		constexpr __m256d zero = { 0.0, 0.0, 0.0, 0.0 };
		const __m256d one = { 1.0, 1.0, 1.0, 1.0 };
		const auto* oneptr = reinterpret_cast<const double*>(&one);
		for (integer d = 0; d != NDIM; ++d) {
			const auto mask1_ = _mm256_cmp_pd(v[d], zero, _CMP_NEQ_UQ);
			const auto mask2_ = _mm256_cmp_pd(v[d], zero, _CMP_EQ_UQ);
			const auto mask1 = *(reinterpret_cast<const __m256i *>(&mask1_));
			const auto mask2 = *(reinterpret_cast<const __m256i *>(&mask2_));
			auto ptr = reinterpret_cast<const double*>(v.data() + d);
			auto a = _mm256_maskload_pd(ptr, mask1);
			const auto* aptr = reinterpret_cast<const double*>(&a);
			const auto b = _mm256_maskload_pd(oneptr, mask2);
			a = _mm256_div_pd(one, _mm256_add_pd(a, b));
			c.v[d] = _mm256_maskload_pd(aptr, mask1);
		}
		//	printf( "%e %e %e %e %e %e %e %e\n", c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);

		return c;
	}
	simd_vector() = default;
	inline ~simd_vector() = default;
	simd_vector(const simd_vector& other) {
		v[0] = other.v[0];
		v[1] = other.v[1];
	}
	inline simd_vector(double d) {
		for (integer i = 0; i != SIMD_SIZE; ++i) {
			v[i] = _mmx_set_pd(d);
		}
	}
	inline double sum() const {
		double r = ZERO;
		for (integer i = 0; i != simd_len; ++i) {
			r += (*this)[i];
		}
		return r;
	}
	inline simd_vector& operator=(const simd_vector& other) {
		v[0] = other.v[0];
		v[1] = other.v[1];
		return *this;
	}
	inline simd_vector operator+(const simd_vector& other) const {
		simd_vector r;
		for (integer i = 0; i != SIMD_SIZE; ++i) {
//_mm_empty();
			r.v[i] = _mmx_add_pd(v[i], other.v[i]);
		}
		return r;
	}
	inline simd_vector operator-(const simd_vector& other) const {
		simd_vector r;
		for (integer i = 0; i != SIMD_SIZE; ++i) {
//_mm_empty();
			r.v[i] = _mmx_sub_pd(v[i], other.v[i]);
		}
		return r;
	}
	inline simd_vector operator*(const simd_vector& other) const {
		simd_vector r;
		for (integer i = 0; i != SIMD_SIZE; ++i) {
//_mm_empty();
			r.v[i] = _mmx_mul_pd(v[i], other.v[i]);
		}
		return r;
	}
	inline simd_vector operator/(const simd_vector& other) const {
		simd_vector r;
		for (integer i = 0; i != SIMD_SIZE; ++i) {
//_mm_empty();
			r.v[i] = _mmx_div_pd(v[i], other.v[i]);
		}
		return r;
	}
	inline simd_vector operator+() const {
		return *this;
	}
	inline simd_vector operator-() const {
		return simd_vector(ZERO) - *this;
	}
	inline simd_vector& operator+=(const simd_vector& other) {
		*this = *this + other;
		return *this;
	}
	inline simd_vector& operator-=(const simd_vector& other) {
		*this = *this - other;
		return *this;
	}
	inline simd_vector& operator*=(const simd_vector& other) {
		*this = *this * other;
		return *this;
	}
	inline simd_vector& operator/=(const simd_vector& other) {
		*this = *this / other;
		return *this;
	}

	inline simd_vector operator*(double d) const {
		const simd_vector other = d;
		return other * *this;
	}
	inline simd_vector operator/(double d) const {
		const simd_vector other = ONE / d;
		return *this * other;
	}

	inline simd_vector operator*=(double d) {
		*this = *this * d;
		return *this;
	}
	inline simd_vector operator/=(double d) {
		*this = *this * (ONE / d);
		return *this;
	}
	inline double& operator[](std::size_t i) {
		double* a = reinterpret_cast<double*>(&v);
		return a[i];
	}
	inline double operator[](std::size_t i) const {
		const double* a = reinterpret_cast<const double*>(&v);
		return a[i];
	}

	double max() const {
		const double a = std::max((*this)[0], (*this)[1]);
		const double b = std::max((*this)[2], (*this)[3]);
		const double c = std::max((*this)[4], (*this)[5]);
		const double d = std::max((*this)[6], (*this)[7]);
		const double e = std::max(a, b);
		const double f = std::max(c, d);
		return std::max(e, f);
	}
	double min() const {
		const double a = std::min((*this)[0], (*this)[1]);
		const double b = std::min((*this)[2], (*this)[3]);
		const double c = std::min((*this)[4], (*this)[5]);
		const double d = std::min((*this)[6], (*this)[7]);
		const double e = std::min(a, b);
		const double f = std::min(c, d);
		return std::min(e, f);
	}
	friend simd_vector sqrt(const simd_vector&);
	friend simd_vector operator*(double, const simd_vector& other);
	friend simd_vector operator/(double, const simd_vector& other);
	friend simd_vector max(const simd_vector& a, const simd_vector& b);
	friend simd_vector min(const simd_vector& a, const simd_vector& b);
	friend simd_vector copysign(const simd_vector& x, const simd_vector& y);
	friend simd_vector blend(const simd_vector& x, const simd_vector& y,
			integer dim);

};

inline simd_vector blend(const simd_vector& x, const simd_vector& y,
		integer dim) {
	simd_vector r;
	constexpr unsigned int mask[2][NDIM] = { { 0xFF00FF00, 0xFFFF0000,
			0x00000000 }, { 0xFF00FF00, 0xFFFF0000, 0xFFFFFFFF } };
	switch (dim) {
	case XDIM:
		r.v[0] = _mm256_blend_pd(x.v[0], y.v[0], 0xA);
		r.v[1] = _mm256_blend_pd(x.v[1], y.v[1], 0xA);
		break;
	case YDIM:
		r.v[0] = _mm256_blend_pd(x.v[0], y.v[0], 0xC);
		r.v[1] = _mm256_blend_pd(x.v[1], y.v[1], 0xC);
		break;
	case ZDIM:
		r.v[0] = _mm256_blend_pd(x.v[0], y.v[0], 0x0);
		r.v[1] = _mm256_blend_pd(x.v[1], y.v[1], 0xF);
		break;
	}
	return r;
}

inline simd_vector copysign(const simd_vector& x, const simd_vector& y) {
	static std::uint64_t m1 = 0x8000000000000000;
	static std::uint64_t m2 = ~m1;
	static const double* m1a = reinterpret_cast<double*>(&m1);
	static const double* m2a = reinterpret_cast<double*>(&m2);
	static const simd_vector mask1(*m1a);
	static const simd_vector mask2(*m2a);
	simd_vector z;
	for (integer i = 0; i != SIMD_SIZE; ++i) {
		const auto tmp1 = _mmx_and_pd(y.v[i], mask1.v[i]);
		const auto tmp2 = _mmx_and_pd(x.v[i], mask2.v[i]);
		z.v[i] = _mmx_or_pd(tmp1, tmp2);
	}
	return z;
}

inline simd_vector sqrt(const simd_vector& vec) {
	simd_vector r;
	for (integer i = 0; i != SIMD_SIZE; ++i) {
		//_mm_empty();
		r.v[i] = _mmx_sqrt_pd(vec.v[i]);
	}
	return r;
}

inline simd_vector operator*(double d, const simd_vector& other) {
	const simd_vector a = d;
	return a * other;
}

inline simd_vector operator/(double d, const simd_vector& other) {
	const simd_vector a = d;
	return a / other;
}

inline void simd_pack(simd_vector* dest, double* src, integer src_len,
		integer pos) {
	for (integer i = 0; i != src_len; ++i) {
		dest[i][pos] = src[i];
	}
}

inline void simd_unpack(double* dest, simd_vector* src, integer src_len,
		integer pos) {
	for (integer i = 0; i != src_len; ++i) {
		dest[i] = src[i][pos];
	}
}

inline simd_vector max(const simd_vector& a, const simd_vector& b) {
	simd_vector r;
	for (integer i = 0; i != SIMD_SIZE; ++i) {
		//_mm_empty();
		r.v[i] = _mmx_max_pd(a.v[i], b.v[i]);
	}
	return r;
}

inline simd_vector min(const simd_vector& a, const simd_vector& b) {
	simd_vector r;
	for (integer i = 0; i != SIMD_SIZE; ++i) {
		//_mm_empty();
		r.v[i] = _mmx_min_pd(a.v[i], b.v[i]);
	}
	return r;
}

inline simd_vector abs(const simd_vector& a) {
	return max(a, -a);
}

#else

#include <hpx/parallel/traits/vector_pack_type.hpp>
#include <hpx/runtime/serialization/datapar.hpp>

#if defined(Vc_HAVE_AVX512F)
using simd_vector = Vc::datapar<double, Vc::datapar_abi::avx512>;
using v4sd = Vc::datapar<double, Vc::datapar_abi::avx>;
#elif defined(Vc_HAVE_AVX)
using simd_vector = typename hpx::parallel::traits::vector_pack_type<double, 8>::type;
using v4sd = Vc::datapar<double, Vc::datapar_abi::avx>;
#else
using simd_vector = typename hpx::parallel::traits::vector_pack_type<double, 8>::type;
using v4sd = typename hpx::parallel::traits::vector_pack_type<double, 4>::type;
#endif

#endif

#endif /* SIMD_VECTOR_HPP_ */
