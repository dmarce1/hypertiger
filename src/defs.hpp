/*
 * types.hpp
 *
 *  Created on: May 26, 2015
 *      Author: dmarce1
 */

#include <hpx/config.hpp>


#define icoarse(i,j,k) ((k) * (NX * NX / 4) + (j) * (NX / 2) + (i))

#ifndef _____DEFS_____HPP
#define _____DEFS_____HPP


#define abort_error() printf( "Error in %s on line %i\n", __FILE__, __LINE__); abort()

#define USE_SIMD

#if !defined(OCTOTIGER_FORCEINLINE)
#   if defined(__NVCC__) || defined(__CUDACC__)
#       define OCTOTIGER_FORCEINLINE inline
#   elif defined(_MSC_VER)
#       define OCTOTIGER_FORCEINLINE __forceinline
#   elif defined(__GNUC__)
#       define OCTOTIGER_FORCEINLINE inline __attribute__ ((__always_inline__))
#   else
#       define OCTOTIGER_FORCEINLINE inline
#   endif
#endif

#include "real.hpp"
typedef long long int integer;

typedef unsigned char byte;

#include <array>
#include <iostream>


constexpr integer NDIM = 3;

constexpr integer INX = 8;
constexpr integer BW = 2;
constexpr integer NX = INX + 2 * BW;
constexpr integer N3 = NX * NX * NX;
constexpr integer DN[NDIM] = { 1, NX, NX * NX };

constexpr integer NDIR = 27;

constexpr integer XDIM = 0;
constexpr integer YDIM = 1;
constexpr integer ZDIM = 2;

constexpr integer FXM = 0;
constexpr integer FXP = 1;
constexpr integer FYM = 2;
constexpr integer FYP = 3;
constexpr integer FZM = 4;
constexpr integer FZP = 5;

constexpr integer NFACE = 2 * NDIM;
constexpr integer NVERTEX = 8;
constexpr integer NCHILD = 8;

constexpr real ZERO = real(0);
constexpr real ONE = real(1);
constexpr real TWO = real(2);
constexpr real THREE = real(3);
constexpr real FOUR = real(4);

constexpr real HALF = real(real(1) / real(2));
constexpr real SIXTH = real(real(1) / real(6));
constexpr real TWELFTH = real(real(1) / real(12));

constexpr real cfl = real(real(2) / real(15));
constexpr integer NRK = 2;
constexpr real rk_beta[2] = { ONE, HALF };

template<typename T>
constexpr inline T sqr(T const& val) {
	return val * val;
}

template<typename T>
constexpr inline T cube(T const& val) {
	return val * val * val;
}

template<typename T>
constexpr inline T average(T const& s1, T const& s2) {
	return 0.5 * (s1 + s2);
}
;

template<typename T>
inline void inplace_average(T& s1, T& s2) {
	s1 = s2 = average(s1, s2);
}
;

template<typename T>
std::size_t write(std::ostream& strm, T && t) {
	typedef typename std::decay<T>::type output_type;
	strm.write(reinterpret_cast<char const*>(&t), sizeof(output_type));
	return sizeof(output_type);
}

template<typename T>
std::size_t write(std::ostream& strm, T* t, std::size_t size) {
	strm.write(reinterpret_cast<char const*>(t), sizeof(T) * size);
	return sizeof(T) * size;
}

template<typename T>
std::size_t read(std::istream& strm, T & t) {
	typedef typename std::decay<T>::type input_type;
	strm.read(reinterpret_cast<char*>(&t), sizeof(input_type));
	return sizeof(input_type);
}

template<typename T>
std::size_t read(std::istream& strm, T* t, std::size_t size) {
	strm.read(reinterpret_cast<char*>(t), sizeof(T) * size);
	return sizeof(T) * size;
}

#endif /* TYPES_HPP_ */
